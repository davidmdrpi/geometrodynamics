"""Separately report the post-freeze scalar-support uniqueness extension."""

import argparse
import json
from pathlib import Path

from geometrodynamics.waves import scalar_esu_uniqueness as su


def failure_controls():
    good={key:True for key in su.REQUIRED_CHECKS}
    rows=[]
    for key in su.REQUIRED_CHECKS:
        for mode in ("missing","failed"):
            checks=good.copy()
            if mode=="missing":del checks[key]
            else:checks[key]=False
            result=su.verdict(checks)
            rows.append(dict(key=key,mode=mode,passed=result["uniqueness"]=="UNRESOLVED"
                             and key in result["failed_checks"]))
    return rows


def run_probe():
    exact=su.exact_certificate()
    pointwise=su.independent_pointwise_checks()
    rank=su.rank_controls()
    failures=failure_controls()
    checks=dict(
        coefficient_expansion=exact["residuals"]["three_conditions"]=="0"
            and all(r["wave_expansion_absolute"]/r["wave_scale"]<1e-10 for r in pointwise),
        rank_and_alignment=exact["all_zero"] and rank["max_rank"]==3
            and rank["alignment_max_relative"]<1e-12,
        regularity_contradiction=exact["residuals"]["regularity_contradiction"]=="0",
        independent_momentum_route=exact["residuals"]["momentum_reciprocal"]=="0"
            and exact["residuals"]["constant_B_wave_residual"]=="0"
            and all(r["momentum_absolute"]<1e-10 for r in pointwise),
        degenerate_controls=rank["zero_B_nonzero_Bdd_matrix_norm"]==0
            and rank["transverse_velocity_matrix_norm"]>.1
            and rank["transverse_acceleration_matrix_norm"]>.1
            and rank["trace_matched_identity_residual"]>.1,
        failure_paths=all(row["passed"] for row in failures)
            and su.verdict({"unrelated":True})["uniqueness"]=="UNRESOLVED",
    )
    checks={key:bool(value) for key,value in checks.items()}
    return dict(post_freeze=True,seed=su.SEED,
                original_freeze="de55f3f3175adafcf2cb760e0aef5767ca3e5016",
                exact=exact,pointwise=pointwise,rank_controls=rank,failure_controls=failures,
                checks=checks,checks_passed=not su.failed_checks(checks),verdict=su.verdict(checks),
                original_freeze_or_archive_changed=False,
                review_clarification="Alignment requires B!=0; B=0 forces Bd=0 but leaves Bdd unrestricted pointwise.",
                proof_scope="one smooth real conformal scalar; exact pointwise round ESU; no other matter; no parity needed for uniqueness")


def render(report):
    lines=["# Post-freeze scalar-support uniqueness","",
           "This extension was not a prediction in freeze de55f3f.","",
           "**"+report["verdict"]["uniqueness"]+"**","",
           "Alignment requires B != 0. At B = 0 only Bd = 0 follows from the pointwise matrix equation.",
           "The independent zero-momentum route also excludes every nonconstant smooth support.","",
           "| Supplementary gate | Pass |","|---|---|"]
    lines += [f"| {key} | {report['checks'].get(key,False)} |" for key in su.REQUIRED_CHECKS]
    if report["verdict"]["failed_checks"]:
        lines += ["","Failed or missing: "+", ".join(report["verdict"]["failed_checks"])]
    lines += ["","The full component-boundary and zero-slice continuation proofs are in",
              "docs/scalar_esu_uniqueness.md. Numerical controls are not a proof by search.","",
              "The even homogeneous family remains outside the imposed odd sector.",
              "Support selection and Phi selection are NOT_DERIVED; the causality gate remains OPEN.",""]
    return "\n".join(lines)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path)
    args=parser.parse_args(argv)
    report=run_probe()
    report["checks_passed"]=not su.failed_checks(report.get("checks",{}))
    report["verdict"]=su.verdict(report.get("checks",{}))
    summary=render(report)
    if args.output_dir:
        args.output_dir.mkdir(parents=True,exist_ok=True)
        (args.output_dir/"uniqueness.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
        (args.output_dir/"uniqueness.md").write_text(summary)
    print(summary)
    return 0 if report["checks_passed"] else 1


if __name__=="__main__":
    raise SystemExit(main())
