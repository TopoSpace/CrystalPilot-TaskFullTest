"""Export inspectable CIF/FCF files, uncertainties, periodic graphs and plots."""
import csv
import itertools
import json
import math
from collections import Counter, deque
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from framework import (ROOT, CELL, INV, VOLUME, OPS, orbit, load_data,
                       structure_factor, fourier_map, SCATTER)
from refine import Model, metrics
from restraints import Restraints


def symmetry_text(op):
    rows = []
    for row in op:
        text = ""
        for c, letter in zip(row, "xyz"):
            if c:
                text += ("+" if c > 0 else "-") + letter
        rows.append(text.lstrip("+") or "0")
    return ",".join(rows)


def cell_header():
    return [
        "_cell_length_a 39.1900", "_cell_length_b 39.1900",
        "_cell_length_c 16.6100", "_cell_angle_alpha 90",
        "_cell_angle_beta 90", "_cell_angle_gamma 120",
        f"_cell_volume {VOLUME:.4f}",
        "_space_group_name_H-M_alt 'P 6/m m m'",
        "_space_group_IT_number 191",
        "_space_group_crystal_system hexagonal",
        "_diffrn_radiation_wavelength 0.68883",
        "_diffrn_ambient_temperature ?",
        "_cell_measurement_temperature ?",
        "_chemical_formula_sum ?",
        "_chemical_formula_weight ?",
        "_cell_formula_units_Z ?",
        "loop_", "_space_group_symop_id", "_space_group_symop_operation_xyz",
    ] + [f"{i+1} '{symmetry_text(op)}'" for i, op in enumerate(OPS)]


def cif_number(value, sigma=None):
    if sigma is None or sigma < 1e-10 or not np.isfinite(sigma):
        return f"{value:.8f}"
    sigma_digits = 1-int(math.floor(math.log10(sigma)))
    rounded_sigma = round(sigma, sigma_digits)
    # Preserve special-position identities and the Cartesian/CIF tensor
    # transformation at interchange precision, even when uncertainties are large.
    esd = round(rounded_sigma*10**8)
    return f"{value:.8f}({esd})"


def uncertainties(name, atoms, log, data, fc):
    initial = json.loads((ROOT / f"{name}_initial.json").read_text())
    normal = dict(np.load(ROOT / f"{name}_normal.npz"))
    selected = (data["d"] >= log["dmin"]) & (data["d"] <= log.get("dmax", 100))
    model = Model(initial, data["hkl"][selected], data["s2"][selected],
                  aniso=log["aniso"], zr_fp=log.get("zr_fp", False))
    p = normal["parameters"]
    reconstructed = model.output_atoms(p)
    assert max(np.max(np.abs(np.array(a["xyz"])-np.array(b["xyz"])))
               for a, b in zip(reconstructed, atoms)) < 1e-10
    scale = log["all_data"]["scale"]
    ic = scale*fc**2
    pp = (np.maximum(data["i"], 0)+2*ic)/3
    weights = 1/(data["sig"]**2+(log["weight"]*pp)**2)
    rss = np.sum(weights[selected]*(data["i"][selected]-ic[selected])**2)
    dof = int(selected.sum())-len(p)
    gof = np.sqrt(rss/dof)
    information = normal["jac"].T @ normal["jac"]
    cov = np.linalg.pinv(information, rcond=1e-12) * gof**2
    esds = {}
    recip_norm = np.linalg.norm(INV, axis=1)
    trans = np.diag(1/recip_norm) @ INV
    parameter_names = []
    for s, a in zip(model.sites, atoms):
        j, n, nu = s["start"], s["n"], s["nu"]
        xyz_cov = s["q"] @ cov[j:j+n, j:j+n] @ s["q"].T
        ub_cif = trans @ s["ub"] @ trans.T
        uj = np.zeros((6, len(p)))
        for k, (a1, a2) in enumerate([(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]):
            uj[k, j+n:j+n+nu] = ub_cif[:, a1, a2]
        uesd = np.sqrt(np.maximum(0, np.diag(uj @ cov @ uj.T)))
        trace_basis = np.trace(s["ub"], axis1=1, axis2=2)/3
        ueq_esd = np.sqrt(trace_basis @ cov[j+n:j+n+nu, j+n:j+n+nu] @ trace_basis)
        esds[a["label"]] = dict(xyz=np.sqrt(np.maximum(0, np.diag(xyz_cov))).tolist(),
                                u_cif=uesd.tolist(), ueq=float(ueq_esd))
        parameter_names += [f'{a["label"]}_position_dof_{k+1}' for k in range(n)]
        parameter_names += [f'{a["label"]}_U_basis_{k+1}' for k in range(nu)]
    if model.fp_index is not None:
        parameter_names.append("Zr_effective_real_scattering_offset")
    parameter_names.append("log_intensity_scale")
    diag = np.sqrt(np.maximum(0, np.diag(cov)))
    corr = cov/np.maximum(np.outer(diag, diag), 1e-30)
    corr_pairs = [(abs(corr[i, j]), parameter_names[i], parameter_names[j],
                   float(corr[i, j])) for i in range(len(p)) for j in range(i)]
    corr_pairs.sort(reverse=True)
    with (ROOT / f"{name}_parameter_uncertainty.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["parameter", "value", "conditional_standard_uncertainty"])
        writer.writerows(zip(parameter_names, p, diag))
    np.savez_compressed(ROOT / f"{name}_covariance.npz", covariance=cov,
                        parameter_names=np.array(parameter_names))
    geometry = Restraints(model)
    ng = len(geometry.calculate(p)[0])
    nactive = int(np.count_nonzero(model.adp_restraints(p)[0]))
    report = dict(n_reflections=int(selected.sum()), n_parameters=len(p),
                  n_restraint_equations=ng+nactive,
                  n_geometry_adp_restraint_equations=ng,
                  n_active_positive_adp_penalties=nactive, gof=float(gof),
                  wr2=float(np.sqrt(rss/np.sum(weights[selected]*data["i"][selected]**2))),
                  normal_matrix_rank=int(np.linalg.matrix_rank(information)),
                  largest_correlations=[dict(a=a, b=b, correlation=r)
                                        for _, a, b, r in corr_pairs[:12]],
                  uncertainty_note="Conditional least-squares uncertainties include "
                  "restraints and fitted residual scale; cell uncertainties, solvent "
                  "model error and systematic scattering-factor error are not included.")
    if model.fp_index is not None:
        report["effective_zr_fp"] = float(p[model.fp_index])
        report["effective_zr_fp_su"] = float(diag[model.fp_index])
    return esds, report, selected


def expand_atoms(atoms):
    expanded = []
    for a in atoms:
        for k, x in enumerate(orbit(a["xyz"])):
            expanded.append(dict(label=f'{a["label"]}_{k+1}', site=a["label"],
                                 element=a["element"], xyz=x, u=a["u"]))
    return expanded


def graph_analysis(atoms):
    expanded = expand_atoms([a for a in atoms if a["element"] != "H"])
    xyz = np.array([a["xyz"] for a in expanded])
    n = len(xyz)
    shifts = np.array(list(itertools.product((-1, 0, 1), repeat=3)))
    points = (shifts[:, None]+xyz[None]).reshape(-1, 3) @ CELL.T
    tree = cKDTree(points)
    adjacent = [[] for _ in range(n)]
    bond_rows = []
    for i, neighbours in enumerate(tree.query_ball_point(xyz @ CELL.T, 2.65)):
        for index in neighbours:
            j = index % n
            translation = shifts[index//n]
            if i == j and not translation.any():
                continue
            ai, aj = expanded[i], expanded[j]
            d = np.linalg.norm(CELL @ (xyz[j]+translation-xyz[i]))
            types = {ai["element"], aj["element"]}
            bond = ((types == {"Zr", "O"} and 1.7 < d < 2.65) or
                    (types == {"C"} and 1.0 < d < 1.85) or
                    (types == {"C", "O"} and 1.0 < d < 1.65))
            if bond:
                adjacent[i].append((j, translation.copy(), d))
                if i < j:
                    bond_rows.append((i, j, translation.copy(), d))
    def components(allowed):
        visited = {}
        groups = []
        cycles = []
        for start in sorted(allowed):
            if start in visited:
                continue
            cid = len(groups)
            group = []
            visited[start] = (cid, np.zeros(3, dtype=int))
            todo = deque([start])
            while todo:
                i = todo.popleft()
                group.append(i)
                ci = visited[i][1]
                for j, tr, _ in adjacent[i]:
                    if j not in allowed:
                        continue
                    cj = ci+tr
                    if j not in visited:
                        visited[j] = (cid, cj)
                        todo.append(j)
                    else:
                        winding = cj-visited[j][1]
                        if winding.any():
                            cycles.append(winding)
            groups.append(group)
        return visited, groups, cycles
    all_map, all_groups, cycles = components(set(range(n)))
    linker_set = {i for i, a in enumerate(expanded)
                  if a["element"] == "C" or a["site"] in ("O4", "O5")}
    node_set = set(range(n))-linker_set
    node_map, nodes, nc = components(node_set)
    linker_map, linkers, lc = components(linker_set)
    connections = set()
    for i in node_set:
        for j, tr, d in adjacent[i]:
            if j not in linker_set:
                continue
            ni, offset_n = node_map[i]
            li, offset_l = linker_map[j]
            delta = tr+offset_n-offset_l
            connections.add((ni, li, *delta.tolist()))
    node_degree = Counter(c[0] for c in connections)
    linker_degree = Counter(c[1] for c in connections)
    zr_coord = [len(adjacent[i]) for i, a in enumerate(expanded) if a["element"] == "Zr"]
    bond_summary = {}
    for i, j, tr, d in bond_rows:
        key = "-".join(sorted((expanded[i]["element"], expanded[j]["element"])))
        bond_summary.setdefault(key, []).append(d)
    bond_summary = {k:dict(count=len(v), minimum=min(v), maximum=max(v), mean=float(np.mean(v)))
                    for k, v in bond_summary.items()}
    adp_eigen = {a["label"]:np.linalg.eigvalsh(a.get("u_cart", np.eye(3)*a["u"])).tolist()
                for a in atoms if a["element"] != "H"}
    result = dict(
        unit_cell_non_hydrogen_atoms=n,
        unit_cell_model_contents=dict(Counter(a["element"] for a in expand_atoms(atoms))),
        connected_components=len(all_groups),
        periodic_connectivity_rank=int(np.linalg.matrix_rank(np.array(cycles))),
        zr_coordination_numbers=zr_coord,
        node_count=len(nodes), linker_count=len(linkers),
        node_compositions=[dict(Counter(expanded[i]["element"] for i in g)) for g in nodes],
        linker_compositions=[dict(Counter(expanded[i]["element"] for i in g)) for g in linkers],
        node_linker_degrees=[node_degree[i] for i in range(len(nodes))],
        linker_node_degrees=[linker_degree[i] for i in range(len(linkers))],
        unique_node_linker_edges=len(connections), bond_summary=bond_summary,
        all_adps_positive=all(min(v) > 0 for v in adp_eigen.values()),
        adp_eigenvalues=adp_eigen)
    assert result["periodic_connectivity_rank"] == 3
    assert sorted(result["node_linker_degrees"]) == [8, 8, 8]
    assert sorted(result["linker_node_degrees"]) == [4]*6
    assert zr_coord == [8]*18
    return result, expanded, bond_rows, nodes, node_map, linkers, linker_map


def write_cif(name, filename, atoms, log, esd, report, selected, data, fc):
    scale = log["all_data"]["scale"]
    m = metrics({k:v[selected] for k,v in data.items()}, fc[selected], scale)
    lines = ["data_NU1000_framework", "# Original independent calculation; PROVISIONAL.",
             "# Not a deposition-ready or externally validated structure."] + cell_header()
    lines += [
        "_chemical_name_common 'NU-1000 framework, provisional'",
        "_computing_structure_solution 'Original Patterson and heavy-atom Fourier algorithms'",
        "_computing_structure_refinement 'Original Python/NumPy/SciPy intensity least squares'",
        "_refine_special_details",
        ";",
        "No crystallographic software or crystallographic library was used.",
        "Neutral-atom Gaussian factors are explicit inputs in framework.py.",
        "Zr real scattering offset is empirical when indicated; f-double-prime=0.",
        "No solvent molecule, guest atom or pore-density correction is included.",
        "All atomic occupancies are fixed to unity; node H atoms are unassigned.",
        "Aromatic H atoms are riding, C-H=0.95 A, Uiso(H)=1.2 Ueq(C).",
        "Aromatic distance/angle/planarity and displacement restraints are applied.",
        "Listed uncertainties are conditional, not complete experimental errors.",
        "The chemical formula and Z are unknown for the solvated measured crystal.",
        f"Model cell contents: C264 H132 O96 Zr18; node hydrogen is incomplete.",
        f"Selected measured d range: {data['d'][selected].min():.6f} to "
        f"{data['d'][selected].max():.6f} A.",
        "See README_zh.md and final_validation.json for limitations and full tests.",
        ";",
        "_refine_ls_structure_factor_coef Fsqd",
        f"_refine_ls_number_reflns {report['n_reflections']}",
        f"_refine_ls_number_parameters {report['n_parameters']}",
        f"_refine_ls_number_restraints {report['n_restraint_equations']}",
        "_reflns_threshold_expression 'I > 2 sigma(I)'",
        f"_reflns_number_gt {m['observed']}",
        f"_refine_ls_R_factor_gt {m['r1_observed']:.6f}",
        f"_refine_ls_R_factor_all {m['r1_all']:.6f}",
        f"_refine_ls_wR_factor_ref {report['wr2']:.6f}",
        f"_refine_ls_goodness_of_fit_ref {report['gof']:.6f}",
        "_refine_ls_weighting_scheme calc",
        "_refine_ls_weighting_details",
        f"'w=1/[sigma(I)^2+({log['weight']}P)^2], P=[max(I,0)+2Ic]/3'",
        f"_reflns_d_resolution_high {data['d'][selected].min():.8f}",
        f"_reflns_d_resolution_low {data['d'][selected].max():.8f}",
        "loop_", "_atom_type_symbol", "_atom_type_scat_dispersion_real",
        "_atom_type_scat_dispersion_imag",
    ]
    for element in ("C", "H", "O", "Zr"):
        fp = next((a.get("fp", 0) for a in atoms if a["element"] == element), 0)
        lines.append(f"{element} {fp:.8f} 0")
    lines += [
        "loop_", "_atom_site_label", "_atom_site_type_symbol",
        "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z",
        "_atom_site_U_iso_or_equiv", "_atom_site_adp_type", "_atom_site_occupancy",
        "_atom_site_symmetry_multiplicity", "_atom_site_calc_flag",
    ]
    for a in atoms:
        e = esd.get(a["label"], {})
        position = " ".join(cif_number(v, s) for v, s in zip(a["xyz"], e.get("xyz", [None]*3)))
        lines.append(f'{a["label"]} {a["element"]} {position} '
                     f'{cif_number(a["u"], e.get("ueq"))} '
                     f'{"Uani" if "u_cart" in a else "Uiso"} 1.0 '
                     f'{len(orbit(a["xyz"]))} {"calc" if a.get("calculated") else "d"}')
    lines += [
        "loop_", "_atom_site_aniso_label", "_atom_site_aniso_U_11",
        "_atom_site_aniso_U_22", "_atom_site_aniso_U_33", "_atom_site_aniso_U_23",
        "_atom_site_aniso_U_13", "_atom_site_aniso_U_12",
    ]
    trans = np.diag(1/np.linalg.norm(INV, axis=1)) @ INV
    for a in atoms:
        if "u_cart" not in a:
            continue
        umat = trans @ a["u_cart"] @ trans.T
        vals = [umat[i, j] for i, j in [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]]
        values = " ".join(cif_number(v, s) for v, s in zip(vals, esd[a["label"]]["u_cif"]))
        lines.append(a["label"]+" "+values)
    (ROOT / filename).write_text("\n".join(lines)+"\n", encoding="ascii")
    return m


def figures(atoms, graph, fc, data, scale):
    _, expanded, bonds, nodes, node_map, linkers, linker_map = graph
    colors = {"C":"#45484b", "O":"#df443f", "Zr":"#009e92", "H":"#cccccc"}
    radii = {"C":7, "O":11, "Zr":36}
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), constrained_layout=True)
    ax = axes[0]
    # Plot cell-crossing segments in both adjoining images, clipped to the axes.
    segments = []
    for i, j, tr, d in bonds:
        a = expanded[i]["xyz"]
        b = expanded[j]["xyz"]+tr
        for shift in (np.zeros(3), -tr):
            segments.append(np.array([CELL@(a+shift), CELL@(b+shift)])[:, :2])
    ax.add_collection(LineCollection(segments, colors="#a5a7aa", linewidths=.5, alpha=.75))
    xyz = np.array([a["xyz"] for a in expanded]) @ CELL.T
    for element in ("C", "O", "Zr"):
        ids = [i for i, a in enumerate(expanded) if a["element"] == element]
        ax.scatter(xyz[ids, 0], xyz[ids, 1], s=radii[element],
                   color=colors[element], edgecolors="none", label=element, zorder=3)
    corners = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 0]]) @ CELL.T
    ax.plot(corners[:, 0], corners[:, 1], color="#202326", lw=.8)
    ax.set_aspect("equal")
    ax.set_xlim(-21, 41)
    ax.set_ylim(-1, 36)
    ax.set_title("Framework along c (no guests)")
    ax.set_xlabel("Cartesian x / A")
    ax.set_ylabel("Cartesian y / A")
    ax.legend(frameon=False, loc="lower right")
    ax = axes[1]
    group = linkers[0]
    coords = np.array([CELL@(expanded[i]["xyz"]+linker_map[i][1]) for i in group])
    centred = coords-coords.mean(axis=0)
    _, _, vh = np.linalg.svd(centred)
    projected = centred @ vh[:2].T
    local = {j:i for i, j in enumerate(group)}
    lines = []
    for i, j, tr, d in bonds:
        if i in local and j in local:
            lines.append([projected[local[i]], projected[local[j]]])
    ax.add_collection(LineCollection(lines, colors="#666a6d", linewidths=1.3))
    for element in ("C", "O"):
        ids = [k for k, i in enumerate(group) if expanded[i]["element"] == element]
        ax.scatter(projected[ids, 0], projected[ids, 1],
                   color=colors[element], s=35 if element == "C" else 65,
                   edgecolors="white", linewidths=.4, zorder=3)
    ax.set_aspect("equal")
    ax.autoscale()
    ax.margins(.12)
    ax.set_title("Recovered C44O8 linker (principal-plane view)")
    ax.set_xlabel("Principal coordinate / A")
    ax.set_ylabel("Principal coordinate / A")
    fig.savefig(ROOT / "framework_views.png", dpi=190)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
    obs = data["i"] > 2*data["sig"]
    fo = np.sqrt(data["i"][obs])
    fca = np.sqrt(scale)*np.abs(fc[obs])
    sc = axes[0].scatter(fo, fca, c=np.minimum(data["d"][obs], 5), s=9,
                         cmap="viridis", alpha=.65)
    limit = max(fo.max(), fca.max())*1.02
    axes[0].plot([0, limit], [0, limit], color="#888888", lw=1)
    axes[0].set(xlabel="Observed |F| (I > 2 sigma)", ylabel="Calculated |F|",
                title="Unmasked all-data model", xlim=(0, limit), ylim=(0, limit))
    fig.colorbar(sc, ax=axes[0], label="d / A (clipped at 5)")
    audit = json.loads((ROOT / "resolution_audit.json").read_text())
    rows = [r for r in audit["shells"] if r["cc_half"] is not None]
    x = [(r["d_min"]+r["d_max"])/2 for r in rows]
    axes[1].plot(x, [r["cc_half"] for r in rows], "o-", color="#008f83")
    axes[1].axhline(0, color="#888888", lw=.8)
    axes[1].set(xlabel="Resolution shell midpoint / A", ylabel="CC1/2",
                title="Independent half-dataset correlation", xlim=(.65, 4.2),
                ylim=(-.2, 1.05))
    fig.savefig(ROOT / "data_quality.png", dpi=190)
    plt.close(fig)


def run():
    data = load_data()
    reports = {}
    exports = [("final_all_data", "nu1000_framework.cif"),
               ("highangle_check", "nu1000_highangle_check.cif"),
               ("no_dispersion", "nu1000_no_dispersion_check.cif")]
    all_atoms = {}
    for name, filename in exports:
        atoms = json.loads((ROOT / f"{name}_atoms.json").read_text())
        log = json.loads((ROOT / f"{name}_log.json").read_text())
        fc = structure_factor(data["hkl"], atoms)
        esd, report, selected = uncertainties(name, atoms, log, data, fc)
        report["refined_subset"] = write_cif(name, filename, atoms, log, esd,
                                            report, selected, data, fc)
        report["all_input_data"] = log["all_data"]
        graph = graph_analysis(atoms)
        report["geometry"] = graph[0]
        reports[name] = report
        all_atoms[name] = atoms
        if name != "final_all_data":
            continue
        scale = log["all_data"]["scale"]
        intensity_calc = scale*fc**2
        with (ROOT / "nu1000_calculated.fcf").open("w", encoding="ascii") as f:
            f.write("data_NU1000_calculated\n")
            f.write("\n".join(cell_header())+"\n")
            f.write("\nloop_\n_refln_index_h\n_refln_index_k\n_refln_index_l\n"
                    "_refln_F_squared_meas\n_refln_F_squared_sigma\n"
                    "_refln_F_squared_calc\n_refln_phase_calc\n")
            for h, io, sig, ic, ph in zip(data["hkl"], data["i"], data["sig"],
                                        intensity_calc, np.where(fc >= 0, 0, 180)):
                f.write(f"{h[0]} {h[1]} {h[2]} {io:.9g} {sig:.9g} {ic:.9g} {ph}\n")
        with (ROOT / "final_reflections.csv").open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["h", "k", "l", "d_A", "I_observed", "sigma_I",
                             "I_calculated", "Fcal_electron_units",
                             "multiplicity", "observed_gt_2sigma", "used_in_refinement"])
            writer.writerows([*h.tolist(), dd, io, si, ic, ff, int(mult),
                              bool(io > 2*si), True] for h, dd, io, si, ic, ff, mult in
                             zip(data["hkl"], data["d"], data["i"], data["sig"],
                                 intensity_calc, fc, data["multiplicity"]))
        expanded = expand_atoms(atoms)
        with (ROOT / "nu1000_unitcell.xyz").open("w", encoding="ascii") as f:
            f.write(str(len(expanded))+"\n")
            lattice = " ".join(f"{v:.8f}" for v in CELL.T.ravel())
            f.write(f'Lattice="{lattice}" Properties=species:S:1:pos:R:3 '
                    'pbc="T T T" provisional_framework_no_guests=true\n')
            for a in expanded:
                cart = CELL @ a["xyz"]
                f.write(f'{a["element"]} {cart[0]:.8f} {cart[1]:.8f} {cart[2]:.8f}\n')
        fo = np.sqrt(np.maximum(data["i"], 0)/scale)
        use = data["d"] >= 1.2
        maps = {}
        for kind, coef in [("fo", fo*np.sign(fc)), ("difference", fo*np.sign(fc)-fc),
                           ("two_fo_minus_fc", 2*fo*np.sign(fc)-fc)]:
            density = fourier_map(data["hkl"][use], coef[use], shape=(128, 128, 64))
            maps[kind] = density.astype(np.float32)
            report[kind+"_map"] = dict(dmin=1.2, maximum=float(density.max()),
                                      minimum=float(density.min()),
                                      rms=float(np.sqrt(np.mean(density**2))))
        np.savez_compressed(ROOT / "final_maps.npz", **maps, cell=CELL, units="electron/A^3",
                            note="Unmasked Fourier sums; F000 and unmeasured reflections set to zero")
        figures(atoms, graph, fc, data, scale)
    sensitivity = {}
    main = all_atoms["final_all_data"]
    for name in ("highangle_check", "no_dispersion"):
        changes = {}
        for a, b in zip(main, all_atoms[name]):
            if a["element"] == "H":
                continue
            diff = np.array(a["xyz"])-b["xyz"]
            diff -= np.round(diff)
            changes[a["label"]] = float(np.linalg.norm(CELL @ diff))
        sensitivity[name] = dict(displacements_A=changes, maximum_A=max(changes.values()),
                                  rms_A=float(np.sqrt(np.mean(np.array(list(changes.values()))**2))))
    reports["sensitivity"] = sensitivity
    reports["holdout_refinement"] = json.loads((ROOT / "holdout_check_log.json").read_text())
    reports["status"] = "Provisional solved framework; publication-grade refinement NOT achieved."
    (ROOT / "final_validation.json").write_text(json.dumps(reports, indent=2))
    (ROOT / "scattering_constants.json").write_text(json.dumps(SCATTER, indent=2))
    print(json.dumps({k:{key:value for key, value in v.items() if key != "geometry"}
                      if k in [a for a, b in exports] else v
                      for k, v in reports.items()}, indent=2))


if __name__ == "__main__":
    run()
