"""Limited, self-written CIF syntax and independent structure-factor checks."""
from pathlib import Path
import ast
import collections
import hashlib
import json
import re
import shlex
import sys
import numpy as np
import pandas as pd
from refine import SCATTER

OUT=Path(__file__).resolve().parent


def cif_tokens(path):
    lines=Path(path).read_text(encoding="ascii").splitlines()
    tokens=[]
    i=0
    while i<len(lines):
        line=lines[i]
        if line.startswith(";"):
            block=[line[1:]]
            i+=1
            while i<len(lines) and not lines[i].startswith(";"):
                block.append(lines[i])
                i+=1
            assert i<len(lines),"Unterminated semicolon text"
            tokens.append("\n".join(block))
        else:
            tokens.extend(shlex.split(line,comments=True,posix=True))
        i+=1
    return tokens


def parse_cif(path):
    tokens=cif_tokens(path)
    scalars={}
    loops=[]
    i=0
    while i<len(tokens):
        token=tokens[i]
        if token.lower().startswith("data_"):
            i+=1
        elif token.lower()=="loop_":
            i+=1
            names=[]
            while i<len(tokens) and tokens[i].startswith("_"):
                names.append(tokens[i])
                i+=1
            values=[]
            while i<len(tokens) and not (tokens[i].startswith("_") or
                    tokens[i].lower()=="loop_" or tokens[i].lower().startswith("data_")):
                values.append(tokens[i])
                i+=1
            assert names and len(values)%len(names)==0,(names,len(values))
            loops.append(pd.DataFrame(np.array(values).reshape(-1,len(names)),columns=names))
        elif token.startswith("_"):
            assert token not in scalars
            scalars[token]=tokens[i+1]
            i+=2
        else:
            raise ValueError(f"Unexpected CIF token: {token}")
    return scalars,loops


def num(s):
    return float(re.sub(r"\(\d+\)$","",str(s)))


def operation(text):
    signs=[]
    shifts=[]
    for axis,part in zip("xyz",text.split(",")):
        assert re.fullmatch(r"-?[xyz](\+1/2)?",part)
        signs.append(-1 if part.startswith("-") else 1)
        shifts.append(.5 if "+1/2" in part else 0.)
        assert axis in part
    return np.array(signs),np.array(shifts)


def check():
    scalars,loops=parse_cif(OUT/"alanine.cif")
    _,floops=parse_cif(OUT/"alanine_reflections.fcf")
    atoms=next(t for t in loops if "_atom_site_label" in t)
    aniso=next(t for t in loops if "_atom_site_aniso_label" in t)
    symmetry=next(t for t in loops if "_space_group_symop_id" in t)
    reflections=floops[0]
    assert len(atoms)==13 and len(aniso)==6 and len(symmetry)==4
    assert len(reflections)==int(scalars["_refine_ls_number_reflns"])==623
    assert collections.Counter(atoms["_atom_site_type_symbol"])=={"C":3,"H":7,"N":1,"O":2}
    assert all(atoms["_atom_site_occupancy"].map(num)==1)
    cell=np.array([num(scalars[f"_cell_length_{a}"]) for a in "abc"])
    h=reflections[["_refln_index_h","_refln_index_k","_refln_index_l"]].astype(int).values
    q=h/cell
    s2=np.sum(q*q,axis=1)/4
    fc=np.zeros(len(h),complex)
    umin=[]
    for _,row in atoms.iterrows():
        xyz=np.array([num(row[f"_atom_site_fract_{a}"]) for a in "xyz"])
        if row["_atom_site_adp_type"]=="Uani":
            ur=aniso[aniso["_atom_site_aniso_label"]==row["_atom_site_label"]].iloc[0]
            u=np.array([[num(ur["_atom_site_aniso_U_11"]),num(ur["_atom_site_aniso_U_12"]),num(ur["_atom_site_aniso_U_13"])],
                        [num(ur["_atom_site_aniso_U_12"]),num(ur["_atom_site_aniso_U_22"]),num(ur["_atom_site_aniso_U_23"])],
                        [num(ur["_atom_site_aniso_U_13"]),num(ur["_atom_site_aniso_U_23"]),num(ur["_atom_site_aniso_U_33"])]])
            umin.append(np.linalg.eigvalsh(u).min())
        else:
            u=np.eye(3)*num(row["_atom_site_U_iso_or_equiv"])
        a,b,c=SCATTER[row["_atom_site_type_symbol"]]
        scattering=np.sum([aa*np.exp(-bb*s2) for aa,bb in zip(a,b)],axis=0)+c
        for text in symmetry["_space_group_symop_operation_xyz"]:
            sg,tr=operation(text)
            xx=xyz*sg+tr
            uu=u*sg[:,None]*sg[None,:]
            thermal=np.exp(-2*np.pi**2*np.einsum("hi,ij,hj->h",q,uu,q))
            fc+=scattering*thermal*np.exp(2j*np.pi*(h@xx))
    assert min(umin)>0
    fc2=abs(fc)**2
    saved=reflections["_refln_F_squared_calc"].astype(float).values
    factor=fc2@saved/(fc2@fc2)
    fracerr=np.linalg.norm(factor*fc2-saved)/np.linalg.norm(saved)
    assert fracerr<.001
    io=reflections["_refln_F_squared_meas"].astype(float).values
    sig=reflections["_refln_F_squared_sigma"].astype(float).values
    fo=np.sqrt(np.maximum(io,0))
    calc=np.sqrt(factor*fc2)
    r1=sum(abs(fo-calc))/sum(fo)
    robs=sum(abs(fo[io>2*sig]-calc[io>2*sig]))/sum(fo[io>2*sig])
    assert abs(r1-num(scalars["_refine_ls_R_factor_all"]))<1e-4
    assert abs(robs-num(scalars["_refine_ls_R_factor_gt"]))<1e-4
    hkl=np.loadtxt(OUT/"alanine_merged.hkl")
    assert np.all(hkl[-1]==0) and len(hkl)-1==len(h)
    assert np.array_equal(hkl[:-1,:3],h)
    assert max(abs(hkl[:-1,3]-io))<=.005001
    allowed={"numpy","scipy","matplotlib","pandas","PIL"}
    local={p.stem for p in OUT.glob("*.py")}
    imports=set()
    for p in OUT.glob("*.py"):
        for node in ast.walk(ast.parse(p.read_text())):
            if isinstance(node,ast.Import):
                imports.update(n.name.split(".")[0] for n in node.names)
            elif isinstance(node,ast.ImportFrom) and node.module:
                imports.add(node.module.split(".")[0])
    assert imports <= allowed|local|set(sys.stdlib_module_names),imports-allowed-local-set(sys.stdlib_module_names)
    report={"CIF_syntax_subset_checks":"passed","atoms":13,"anisotropic_atoms":6,
            "reflections":len(h),"independent_Fc2_relative_norm_error_from_rounded_CIF":float(fracerr),
            "R1_all_recomputed_from_CIF":float(r1),"R1_observed_recomputed_from_CIF":float(robs),
            "only_allowed_or_standard_library_imports":True,
            "not_performed":"External CIF dictionary validation or checkCIF/PLATON."}
    (OUT/"output_checks.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report),flush=True)
    return report


def provenance():
    inputs=["environment.json","inputs/alanine/pgw240033_Mo.par",
            "inputs/alanine/original/pgw240033_Mo_cracker.par",
            "inputs/alanine/pgw240033_Mo.run",
            "inputs/alanine/RED_online_log.txt",
            "inputs/alanine/frames/jpg/pgw240033_Mo_diffimg_1.jpg",
            "inputs/alanine/HyPixArc100_290622_calib_170723.modulegeo",
            "inputs/alanine/pgw240033_Mo.CAP_shape",
            "inputs/alanine/expinfo/pgw240033_Mo_crystal.ini",
            "inputs/alanine/expinfo/pgw240033_Mo_sample.ini",
            "inputs/alanine/expinfo/pgw240033_Mo_datacoll.ini",
            "inputs/alanine/pgw240033_Mo_ccd.sum"]
    records=[]
    for name in inputs:
        p=OUT.parent/name
        records.append({"file":name,"sha256":hashlib.sha256(p.read_bytes()).hexdigest()})
    records.extend({"file":f"inputs/alanine/frames/pgw240033_Mo_{r['run']}_{r['frame']}.rodhypix",
                    "sha256":r["sha256"]} for r in json.loads((OUT/"frame_manifest.json").read_text()))
    (OUT/"input_provenance.json").write_text(json.dumps(records,indent=2))
    lines=[]
    extensions={".py",".json",".csv",".cif",".fcf",".hkl",".xyz",".md",".png"}
    for p in sorted(OUT.iterdir()):
        if p.is_file() and p.suffix in extensions:
            lines.append(hashlib.sha256(p.read_bytes()).hexdigest()+"  "+p.name)
    (OUT/"checksums.sha256").write_text("\n".join(lines)+"\n")


if __name__=="__main__":
    check()
    provenance()
