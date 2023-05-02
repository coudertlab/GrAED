// Copyright 2018 Global Phasing Ltd.
//
// Topo(logy) - restraints (from a monomer library) applied to a model.

#ifndef GEMMI_TOPO_HPP_
#define GEMMI_TOPO_HPP_

#include <map>           // for multimap
#include <ostream>       // for ostream
#include <memory>        // for unique_ptr
#include "chemcomp.hpp"  // for ChemComp
#include "monlib.hpp"    // for MonLib
#include "model.hpp"     // for Residue, Atom
#include "calculate.hpp" // for calculate_angle, calculate_dihedral
#include "polyheur.hpp"  // for are_connected

namespace gemmi {

struct Topo {
  // We have internal pointers in this class (pointers setup in
  // apply_restraints() that point to ResInfo::chemcomp.rt),
  // disable copying this class.
  Topo() = default;
  Topo(Topo const&) = delete;
  Topo& operator=(Topo const&) = delete;

  struct Bond {
    const Restraints::Bond* restr;
    std::array<Atom*, 2> atoms;
    double calculate() const { return atoms[0]->pos.dist(atoms[1]->pos); }
    double calculate_z() const {
      return std::abs(calculate() - restr->value) / restr->esd;
    }
  };
  struct Angle {
    const Restraints::Angle* restr;
    std::array<Atom*, 3> atoms;
    double calculate() const {
      return calculate_angle(atoms[0]->pos, atoms[1]->pos, atoms[2]->pos);
    }
    double calculate_z() const { return angle_z(calculate(), *restr); }
  };
  struct Torsion {
    const Restraints::Torsion* restr;
    std::array<Atom*, 4> atoms;
    double calculate() const {
      return calculate_dihedral(atoms[0]->pos, atoms[1]->pos,
                                atoms[2]->pos, atoms[3]->pos);
    }
    double calculate_z() const {
      return angle_z(calculate(), *restr, 360. / std::max(1, restr->period));
    }
  };
  struct Chirality {
    const Restraints::Chirality* restr;
    std::array<Atom*, 4> atoms;
    double calculate() const {
      return calculate_chiral_volume(atoms[0]->pos, atoms[1]->pos,
                                     atoms[2]->pos, atoms[3]->pos);
    }
    double calculate_z(double ideal_abs_vol, double esd) const {
      double calc = calculate();
      if (restr->sign == ChiralityType::Negative ||
          (restr->sign == ChiralityType::Both && calc < 0))
        ideal_abs_vol *= -1;
      return std::abs(calc - ideal_abs_vol) / esd;
    }
    bool check() const { return !restr->is_wrong(calculate()); }
  };
  struct Plane {
    const Restraints::Plane* restr;
    std::vector<Atom*> atoms;
    bool has(const Atom* atom) const {
      return in_vector(const_cast<Atom*>(atom), atoms);
    }
  };

  enum class RKind { Bond, Angle, Torsion, Chirality, Plane };
  struct Rule {
    RKind rkind;
    size_t index; // index in the respective vector (bonds, ...) in Topo
  };

  struct Link {
    std::string link_id;
    Residue* res1 = nullptr;
    Residue* res2 = nullptr;
    std::vector<Rule> link_rules;
    // altloc and asu are used only for ChainInfo::extras, not for ResInfo::prev
    char alt1 = '\0';
    char alt2 = '\0';
    Asu asu = Asu::Any;
    ChemComp::Group aliasing1 = ChemComp::Group::Null;
    ChemComp::Group aliasing2 = ChemComp::Group::Null;

    // only for polymer links, res1 and res2 must be in the same vector (Chain)
    std::ptrdiff_t res_distance() const { return res1 - res2; }
  };

  struct ResInfo {
    Residue* res;
    // in case of microheterogeneity we may have 2+ previous residues
    std::vector<Link> prev;
    std::vector<std::pair<std::string, ChemComp::Group>> mods;
    ChemComp chemcomp;
    std::vector<Rule> monomer_rules;

    ResInfo(Residue* r) : res(r) {}
    void add_mod(const std::string& m, ChemComp::Group aliasing) {
      if (!m.empty())
        mods.emplace_back(m, aliasing);
    }
  };

  // corresponds to a sub-chain
  struct ChainInfo {
    const Chain& chain_ref;
    std::string subchain_name;
    std::string entity_id;
    bool polymer;
    PolymerType polymer_type;
    std::vector<ResInfo> res_infos;

    ChainInfo(ResidueSpan& subchain, const Chain& chain, const Entity* ent);
    void setup_polymer_links();
    struct RGroup {
      std::vector<ResInfo>::iterator begin, end;
    };
    RGroup group_from(std::vector<ResInfo>::iterator b) const {
      auto e = b + 1;
      while (e != res_infos.end() && e->res->group_key() == b->res->group_key())
        ++e;
      return RGroup{b, e};
    }
  private:
    Link make_polymer_link(const ResInfo& ri1, const ResInfo& ri2) const;
  };

  template<typename T>
  static int has_atom(const Atom* a, const T& t) {
    for (int i = 0; (size_t) i != t.atoms.size(); ++i)
      if (t.atoms[i] == a)
        return i;
    return -1;
  }

  std::ostream* warnings = nullptr;
  std::vector<ChainInfo> chain_infos;
  std::vector<Link> extras;

  // Restraints applied to Model
  std::vector<Bond> bonds;
  std::vector<Angle> angles;
  std::vector<Torsion> torsions;
  std::vector<Chirality> chirs;
  std::vector<Plane> planes;

  std::multimap<const Atom*, Bond*> bond_index;       // indexes both atoms
  std::multimap<const Atom*, Angle*> angle_index;     // only middle atom
  std::multimap<const Atom*, Torsion*> torsion_index; // two middle atoms
  std::multimap<const Atom*, Plane*> plane_index;     // all atoms

  ResInfo* find_resinfo(const Residue* res) {
    for (ChainInfo& ci : chain_infos)
      for (ResInfo& ri : ci.res_infos)
        if (ri.res == res)
          return &ri;
    return nullptr;
  }

  Bond* first_bond_in_link(const Link& link) {
    for (const Rule& rule : link.link_rules)
      if (rule.rkind == RKind::Bond)
        return &bonds[rule.index];
    return nullptr;
  }

  const Restraints::Bond* take_bond(const Atom* a, const Atom* b) const {
    auto range = bond_index.equal_range(a);
    for (auto i = range.first; i != range.second; ++i) {
      const Bond* bond = i->second;
      if ((bond->atoms[0] == b && bond->atoms[1] == a) ||
          (bond->atoms[1] == b && bond->atoms[0] == a))
        return bond->restr;
    }
    return nullptr;
  }

  const Restraints::Angle* take_angle(const Atom* a,
                                      const Atom* b,
                                      const Atom* c) const {
    auto range = angle_index.equal_range(b);
    for (auto i = range.first; i != range.second; ++i) {
      const Angle* ang = i->second;
      if ((ang->atoms[0] == a && ang->atoms[2] == c) ||
          (ang->atoms[0] == c && ang->atoms[2] == a))
        return ang->restr;
    }
    return nullptr;
  }

  const Chirality* get_chirality(const Atom* ctr) const {
    for (const Chirality& chir : chirs)
      if (chir.atoms[0] == ctr)
        return &chir;
    return nullptr;
  }

  double ideal_chiral_abs_volume(const Chirality &ch) const {
    const Restraints::Bond* bond_c1 = take_bond(ch.atoms[0], ch.atoms[1]);
    const Restraints::Bond* bond_c2 = take_bond(ch.atoms[0], ch.atoms[2]);
    const Restraints::Bond* bond_c3 = take_bond(ch.atoms[0], ch.atoms[3]);
    const Restraints::Angle* angle_1c2 = take_angle(ch.atoms[1], ch.atoms[0], ch.atoms[2]);
    const Restraints::Angle* angle_2c3 = take_angle(ch.atoms[2], ch.atoms[0], ch.atoms[3]);
    const Restraints::Angle* angle_3c1 = take_angle(ch.atoms[3], ch.atoms[0], ch.atoms[1]);
    if (bond_c1 && bond_c2 && bond_c3 && angle_1c2 && angle_2c3 && angle_3c1)
      return chiral_abs_volume(bond_c1->value, bond_c2->value, bond_c3->value,
                               angle_1c2->value, angle_2c3->value, angle_3c1->value);
    return std::numeric_limits<double>::quiet_NaN();
  }

  std::vector<Rule> apply_restraints(const Restraints& rt,
                                     Residue& res, Residue* res2,
                                     char altloc='*') {
    std::string altlocs;
    if (altloc == '*') {
      // find all distinct altlocs
      add_distinct_altlocs(res, altlocs);
      if (res2)
        add_distinct_altlocs(*res2, altlocs);
    }
    if (altlocs.empty())
      altlocs += altloc;

    std::vector<Rule> rules;
    for (const Restraints::Bond& bond : rt.bonds)
      for (char alt : altlocs)
        if (Atom* at1 = bond.id1.get_from(res, res2, alt))
          if (Atom* at2 = bond.id2.get_from(res, res2, alt)) {
            rules.push_back({RKind::Bond, bonds.size()});
            bonds.push_back({&bond, {{at1, at2}}});
            if (!at1->altloc && !at2->altloc)
              break;
          }
    for (const Restraints::Angle& angle : rt.angles)
      for (char alt : altlocs)
        if (Atom* at1 = angle.id1.get_from(res, res2, alt))
          if (Atom* at2 = angle.id2.get_from(res, res2, alt))
            if (Atom* at3 = angle.id3.get_from(res, res2, alt)) {
              rules.push_back({RKind::Angle, angles.size()});
              angles.push_back({&angle, {{at1, at2, at3}}});
              if (!at1->altloc && !at2->altloc && !at3->altloc)
                break;
            }
    for (const Restraints::Torsion& tor : rt.torsions)
      for (char alt : altlocs)
        if (Atom* at1 = tor.id1.get_from(res, res2, alt))
          if (Atom* at2 = tor.id2.get_from(res, res2, alt))
            if (Atom* at3 = tor.id3.get_from(res, res2, alt))
              if (Atom* at4 = tor.id4.get_from(res, res2, alt)) {
                rules.push_back({RKind::Torsion, torsions.size()});
                torsions.push_back({&tor, {{at1, at2, at3, at4}}});
                if (!at1->altloc && !at2->altloc &&
                    !at3->altloc && !at4->altloc)
                  break;
          }
    for (const Restraints::Chirality& chir : rt.chirs)
      for (char alt : altlocs)
        if (Atom* at1 = chir.id_ctr.get_from(res, res2, alt))
          if (Atom* at2 = chir.id1.get_from(res, res2, alt))
            if (Atom* at3 = chir.id2.get_from(res, res2, alt))
              if (Atom* at4 = chir.id3.get_from(res, res2, alt)) {
                rules.push_back({RKind::Chirality, chirs.size()});
                chirs.push_back({&chir, {{at1, at2, at3, at4}}});
                if (!at1->altloc && !at2->altloc &&
                    !at3->altloc && !at4->altloc)
                  break;
              }
    for (const Restraints::Plane& plane : rt.planes)
      for (char alt : altlocs) {
        std::vector<Atom*> atoms;
        for (const Restraints::AtomId& id : plane.ids)
          if (Atom* atom = id.get_from(res, res2, alt))
            atoms.push_back(atom);
        if (atoms.size() >= 4) {
          rules.push_back({RKind::Plane, planes.size()});
          planes.push_back({&plane, atoms});
        }
        if (std::all_of(atoms.begin(), atoms.end(),
                        [](Atom* a) { return !a->altloc; }))
          break;
      }
    return rules;
  }

  void apply_restraints_to_residue(ResInfo& ri, const MonLib& monlib) {
    // link restraints
    for (Link& link : ri.prev)
      if (const ChemLink* chem_link = monlib.get_link(link.link_id)) {
        const Restraints* rt = &chem_link->rt;
        // aliases are a new and rarely used thing
        if (link.aliasing1 != ChemComp::Group::Null ||
            link.aliasing2 != ChemComp::Group::Null) {
          std::unique_ptr<Restraints> rt_copy(new Restraints(*rt));
          if (link.aliasing1 != ChemComp::Group::Null) {
            const ChemComp& cc = (&ri + link.res_distance())->chemcomp;
            for (const auto& p : cc.get_aliasing(link.aliasing1).related)
              rt_copy->rename_atom(Restraints::AtomId{1, p.first}, p.second);
          }
          if (link.aliasing2 != ChemComp::Group::Null) {
            const ChemComp& cc = ri.chemcomp;
            for (const auto& p : cc.get_aliasing(link.aliasing2).related)
              rt_copy->rename_atom(Restraints::AtomId{2, p.first}, p.second);
          }
          rt = rt_copy.get();
          rt_storage.push_back(std::move(rt_copy));
        }
        auto rules = apply_restraints(*rt, *link.res1, ri.res);
        vector_move_extend(link.link_rules, std::move(rules));
      }
    // monomer restraints
    auto rules = apply_restraints(ri.chemcomp.rt, *ri.res, nullptr);
    vector_move_extend(ri.monomer_rules, std::move(rules));
  }

  void apply_restraints_to_extra_link(Link& link, const MonLib& monlib) {
    const ChemLink* cl = monlib.get_link(link.link_id);
    if (!cl) {
      err("ignoring link '" + link.link_id + "' as it is not in the monomer library");
      return;
    }
    if (link.alt1 && link.alt2 && link.alt1 != link.alt2)
      err(cat("LINK between different conformers ", link.alt1, " and ", link.alt2, '.'));
    char alt = link.alt1 ? link.alt1 : link.alt2;
    auto rules = apply_restraints(cl->rt, *link.res1, link.res2, alt);
    vector_move_extend(link.link_rules, std::move(rules));
  }

  // Structure is non-const b/c connections may have link_id assigned.
  // Model is non-const b/c we store non-const pointers to residues in Topo.
  // Because of the pointers, don't add or remove residues after this step.
  // Monlib may get modified by addition of extra links from the model.
  void initialize_refmac_topology(Structure& st, Model& model0,
                                  MonLib& monlib, bool ignore_unknown_links=false);

  // This step stores pointers to gemmi::Atom's from model0,
  // so after this step don't add or remove atoms.
  // monlib is needed only for links.
  void finalize_refmac_topology(const MonLib& monlib) {
    for (ChainInfo& chain_info : chain_infos)
      for (ResInfo& ri : chain_info.res_infos)
        apply_restraints_to_residue(ri, monlib);
    for (Link& link : extras)
      apply_restraints_to_extra_link(link, monlib);

    // create indices
    for (Bond& bond : bonds) {
      bond_index.emplace(bond.atoms[0], &bond);
      if (bond.atoms[1] != bond.atoms[0])
        bond_index.emplace(bond.atoms[1], &bond);
    }
    for (Angle& ang : angles)
      angle_index.emplace(ang.atoms[1], &ang);
    for (Torsion& tor : torsions) {
      torsion_index.emplace(tor.atoms[1], &tor);
      if (tor.atoms[1] != tor.atoms[2])
        torsion_index.emplace(tor.atoms[2], &tor);
    }
    for (Plane& plane : planes)
      for (Atom* atom : plane.atoms)
        plane_index.emplace(atom, &plane);
  }

  Link* find_polymer_link(const AtomAddress& a1, const AtomAddress& a2) {
    for (ChainInfo& ci : chain_infos)
      if (a1.chain_name == ci.chain_ref.name && a2.chain_name == ci.chain_ref.name) {
        for (ResInfo& ri : ci.res_infos)
          for (Link& link : ri.prev) {
            assert(link.res1 && link.res2);
            if ((a1.res_id.matches_noseg(*link.res1) &&
                 a2.res_id.matches_noseg(*link.res2)) ||
                (a2.res_id.matches_noseg(*link.res1) &&
                 a1.res_id.matches_noseg(*link.res2)))
              return &link;
          }
      }
    return nullptr;
  }

  GEMMI_COLD void err(const std::string& msg) const {
    if (warnings == nullptr)
      fail(msg);
    *warnings << "Warning: " << msg << std::endl;
  }

private:
  void setup_connection(Connection& conn, Model& model0, MonLib& monlib,
                        bool ignore_unknown_links);
  std::vector<std::unique_ptr<Restraints>> rt_storage;
};

inline Topo::ChainInfo::ChainInfo(ResidueSpan& subchain,
                                  const Chain& chain, const Entity* ent)
  : chain_ref(chain) {
  subchain_name = subchain.at(0).subchain;
  res_infos.reserve(subchain.size());
  if (ent) {
    entity_id = ent->name;
    polymer = ent->entity_type == EntityType::Polymer;
    polymer_type = get_or_check_polymer_type(ent, subchain);
  } else {
    polymer = false;
    polymer_type = PolymerType::Unknown;
  }
  for (Residue& res : subchain)
    res_infos.emplace_back(&res);
}

inline Topo::Link Topo::ChainInfo::make_polymer_link(const Topo::ResInfo& ri1,
                                                     const Topo::ResInfo& ri2) const {
  Link link;
  link.res1 = ri1.res;
  link.res2 = ri2.res;
  assert(&ri1 - &ri2 == link.res_distance());
  if (is_polypeptide(polymer_type)) {
    bool groups_ok = true;
    std::string c = "C";
    std::string n = "N";
    if (!ChemComp::is_peptide_group(ri1.chemcomp.group)) {
      for (const ChemComp::Aliasing& aliasing : ri1.chemcomp.aliases)
        if (ChemComp::is_peptide_group(aliasing.group)) {
          link.aliasing1 = aliasing.group;
          if (const std::string* c_ptr = aliasing.name_from_alias(c))
            c = *c_ptr;
        }
      if (link.aliasing1 == ChemComp::Group::Null)
        groups_ok = false;
    }
    ChemComp::Group n_terminus_group = ri2.chemcomp.group;
    if (!ChemComp::is_peptide_group(ri2.chemcomp.group)) {
      for (const ChemComp::Aliasing& aliasing : ri2.chemcomp.aliases)
        if (ChemComp::is_peptide_group(aliasing.group)) {
          link.aliasing2 = n_terminus_group = aliasing.group;
          if (const std::string* n_ptr = aliasing.name_from_alias(n))
            n = *n_ptr;
        }
      if (link.aliasing2 == ChemComp::Group::Null)
        groups_ok = false;
    }
    const Atom* a1 = ri1.res->find_atom(c, '*', El::C);
    const Atom* a2 = ri2.res->find_atom(n, '*', El::N);
    if (groups_ok && in_peptide_bond_distance(a1, a2)) {
      bool is_cis = ri1.res->is_cis;
      if (n_terminus_group == ChemComp::Group::PPeptide)
        link.link_id = is_cis ? "PCIS" : "PTRANS";
      else if (n_terminus_group == ChemComp::Group::MPeptide)
        link.link_id = is_cis ? "NMCIS" : "NMTRANS";
      else
        link.link_id = is_cis ? "CIS" : "TRANS";
    } else {
      link.link_id = "gap";
    }
  } else if (is_polynucleotide(polymer_type)) {
    std::string o3p = "O3'";
    std::string p = "P";
    bool groups_ok = true;
    if (!ChemComp::is_nucleotide_group(ri1.chemcomp.group)) {
      /* disabled for now
      for (const ChemComp::Aliasing& aliasing : ri1.chemcomp.aliases)
        if (ChemComp::is_nucleotide_group(aliasing.group)) {
          link.aliasing1 = aliasing.group;
          if (const std::string* o3p_ptr = aliasing.name_from_alias(o3p))
            o3p = *o3p_ptr;
        }
      */
      if (link.aliasing1 == ChemComp::Group::Null)
        groups_ok = false;
    }
    if (!ChemComp::is_nucleotide_group(ri2.chemcomp.group)) {
      /* disabled for now
      for (const ChemComp::Aliasing& aliasing : ri2.chemcomp.aliases)
        if (ChemComp::is_nucleotide_group(aliasing.group)) {
          link.aliasing2 = aliasing.group;
          if (const std::string* p_ptr = aliasing.name_from_alias(p))
            p = *p_ptr;
        }
      */
      if (link.aliasing2 == ChemComp::Group::Null)
        groups_ok = false;
    }
    const Atom* a1 = ri1.res->find_atom(o3p, '*', El::O);
    const Atom* a2 = ri2.res->find_atom(p, '*', El::P);
    if (groups_ok && in_nucleotide_bond_distance(a1, a2)) {
      link.link_id = "p";
    } else {
      link.link_id = "gap";
    }
  } else {
    link.link_id = "?";
  }
  return link;
}

inline void Topo::ChainInfo::setup_polymer_links() {
  if (!polymer || res_infos.empty())
    return;
  RGroup prev_group = group_from(res_infos.begin());
  while (prev_group.end != res_infos.end()) {
    RGroup group = group_from(prev_group.end);
    for (auto ri = group.begin; ri != group.end; ++ri)
      for (auto prev_ri = prev_group.begin; prev_ri != prev_group.end; ++prev_ri)
        ri->prev.push_back(make_polymer_link(*prev_ri, *ri));
    prev_group = group;
  }
}

// see comments above the declaration
inline void Topo::initialize_refmac_topology(Structure& st, Model& model0,
                                             MonLib& monlib, bool ignore_unknown_links) {
  // initialize chains and residues
  for (Chain& chain : model0.chains)
    for (ResidueSpan& sub : chain.subchains()) {
      // set Residue::group_idx which is used in Restraints::AtomId::get_from()
      for (size_t i = 0; i != sub.size(); ++i) {
        sub[i].group_idx = 0;
        if (i != 0 && sub[i-1].seqid == sub[i].seqid)
          sub[i].group_idx = sub[i-1].group_idx + 1;
      }
      const Entity* ent = st.get_entity_of(sub);
      chain_infos.emplace_back(sub, chain, ent);
    }
  for (ChainInfo& ci : chain_infos) {
    // copy monomer description
    for (ResInfo& ri : ci.res_infos) {
      auto it = monlib.monomers.find(ri.res->name);
      if (it != monlib.monomers.end())
        ri.chemcomp = it->second;
      else
        err("unknown chemical component " + ri.res->name
            + " in chain " +  ci.chain_ref.name);
    }
    ci.setup_polymer_links();
  }

  // add extra links
  for (Connection& conn : st.connections)
    if (conn.type != Connection::Hydrog) // ignoring hydrogen bonds
      setup_connection(conn, model0, monlib, ignore_unknown_links);

  // Add modifications from standard links. We do it here b/c polymer links
  // could be disabled (link_id = "?") in setup_connection().
  for (ChainInfo& ci : chain_infos)
    for (ResInfo& ri : ci.res_infos)
      for (Link& prev : ri.prev)
        if (const ChemLink* chem_link = monlib.get_link(prev.link_id)) {
          ResInfo* ri_prev = &ri + prev.res_distance();
          ri_prev->add_mod(chem_link->side1.mod, prev.aliasing1);
          ri.add_mod(chem_link->side2.mod, prev.aliasing2);
        }

  for (ChainInfo& chain_info : chain_infos)
    for (ResInfo& ri : chain_info.res_infos) {
      // apply modifications
      for (const auto& modif : ri.mods) {
        if (const ChemMod* chem_mod = monlib.get_mod(modif.first))
          try {
            chem_mod->apply_to(ri.chemcomp, modif.second);
          } catch(std::runtime_error& e) {
            err("failed to apply modification " + chem_mod->id
                + " to " + ri.res->name + ": " + e.what());
          }
        else
          err("modification not found: " + modif.first);
      }
    }
}

// it has side-effects: may modify conn.link_id and add to monlib.links
inline void Topo::setup_connection(Connection& conn, Model& model0, MonLib& monlib,
                                   bool ignore_unknown_links) {
  if (conn.link_id == "gap") {
    Link* polymer_link = find_polymer_link(conn.partner1, conn.partner2);
    if (polymer_link) polymer_link->link_id = "?";  // disable polymer link
    return;
  }

  Link extra;
  CRA cra1 = model0.find_cra(conn.partner1, true);
  CRA cra2 = model0.find_cra(conn.partner2, true);
  if (!cra1.atom || !cra2.atom)
    return;
  extra.res1 = cra1.residue;
  extra.res2 = cra2.residue;
  extra.alt1 = conn.partner1.altloc;
  extra.alt2 = conn.partner2.altloc;
  extra.asu = conn.asu;

  const ChemLink* match = nullptr;

  // If we have link_id find ChemLink by name (and check if it matches).
  if (!conn.link_id.empty()) {
    match = monlib.get_link(conn.link_id);
    if (!match) {
      err("link not found in monomer library: " + conn.link_id);
      return;
    }
    if (match->rt.bonds.empty() ||
        match->rt.bonds[0].id1.atom != conn.partner1.atom_name ||
        match->rt.bonds[0].id2.atom != conn.partner2.atom_name ||
        !monlib.link_side_matches_residue(match->side1, extra.res1->name) ||
        !monlib.link_side_matches_residue(match->side2, extra.res2->name)) {
      err("link from the monomer library does not match: " + conn.link_id);
      return;
    }
  } else {
    // we don't have link_id - use the best matching link (if any)
    auto r = monlib.match_link(*extra.res1, conn.partner1.atom_name,
                               *extra.res2, conn.partner2.atom_name,
                               extra.alt1 ? extra.alt1 : extra.alt2);
    match = r.first;
    if (match && r.second) {
      std::swap(extra.res1, extra.res2);
      std::swap(extra.alt1, extra.alt2);
    }
  }

  // If a polymer link is also given in LINK/struct_conn,
  // use only one of them. If LINK has explicit name (ccp4_link_id),
  // or if it matches residue-specific link from monomer library, use it;
  // otherwise, LINK is repetition of TRANS/CIS, so ignore LINK.
  if (Link* polymer_link = find_polymer_link(conn.partner1, conn.partner2)) {
    if (conn.link_id.empty() && !cif::is_null(polymer_link->link_id) &&
        polymer_link->link_id != "gap" &&
        (!match || (match->side1.comp.empty() && match->side2.comp.empty())))
      return;
    polymer_link->link_id = "?";  // disable polymer link
  }

  if (match) {
    extra.link_id = match->id;
    // add modifications from the link
    find_resinfo(extra.res1)->add_mod(match->side1.mod, extra.aliasing1);
    find_resinfo(extra.res2)->add_mod(match->side2.mod, extra.aliasing2);
  } else {
    if (ignore_unknown_links)
      return;
    // create a new ChemLink and add it to the monomer library
    ChemLink cl;
    cl.side1.comp = extra.res1->name;
    cl.side2.comp = extra.res2->name;
    cl.id = cl.side1.comp + cl.side2.comp;
    cl.name = "auto-" + cl.id;
    bool use_ion = cra1.atom->element.is_metal() || cra2.atom->element.is_metal();
    double ideal_dist = monlib.find_radius(cra1, use_ion) +
                        monlib.find_radius(cra2, use_ion);
    if (!use_ion)
      ideal_dist /= 2;
    cl.rt.bonds.push_back({Restraints::AtomId{1, conn.partner1.atom_name},
                           Restraints::AtomId{2, conn.partner2.atom_name},
                           BondType::Unspec, false,
                           ideal_dist, 0.02,
                           ideal_dist, 0.02});
    monlib.ensure_unique_link_name(cl.id);
    monlib.links.emplace(cl.id, cl);
    extra.link_id = cl.id;
  }
  if (conn.link_id.empty())
    conn.link_id = extra.link_id;
  extras.push_back(extra);
}

} // namespace gemmi
#endif
