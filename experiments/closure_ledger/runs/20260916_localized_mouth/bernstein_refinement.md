# Prospective localized-mouth refinement

Original: **6/8**, unchanged.

Refinement: **6/8**.

```json
{
  "gates": {
    "action": true,
    "momentum": true,
    "hamiltonian": true,
    "seam": true,
    "physical": false,
    "localization": true,
    "controls": true,
    "evidence": false
  },
  "passed": 6,
  "total": 8,
  "metrics": {
    "cases": [
      {
        "L": 3.5,
        "eta": 0.0,
        "H_residual": 4.81229862625554e-08,
        "M_residual": 0.0,
        "minimum_psi": 0.44774691070344474,
        "relative_refinement": 4.116137430007164e-12
      },
      {
        "L": 3.5,
        "eta": 0.1,
        "H_residual": 5.445076944887717e-08,
        "M_residual": 1.0145994423305452e-14,
        "minimum_psi": 0.44774692957070394,
        "relative_refinement": 4.1172740595253905e-12
      },
      {
        "L": 3.5,
        "eta": 0.3,
        "H_residual": 4.4040660091004824e-08,
        "M_residual": 3.043798496398225e-14,
        "minimum_psi": 0.4477470805083947,
        "relative_refinement": 4.113718357607272e-12
      },
      {
        "L": 4.5,
        "eta": 0.0,
        "H_residual": 2.8065835044888132e-08,
        "M_residual": 0.0,
        "minimum_psi": 0.28275831328527573,
        "relative_refinement": 8.788570070989571e-11
      },
      {
        "L": 4.5,
        "eta": 0.1,
        "H_residual": 3.211739252328982e-08,
        "M_residual": 4.200630276569658e-16,
        "minimum_psi": 0.28275831462025464,
        "relative_refinement": 8.788361596991714e-11
      },
      {
        "L": 4.5,
        "eta": 0.3,
        "H_residual": 3.686066835095758e-08,
        "M_residual": 1.2601890697360075e-15,
        "minimum_psi": 0.2827583253000946,
        "relative_refinement": 8.78767360659962e-11
      },
      {
        "L": 5.5,
        "eta": 0.0,
        "H_residual": 2.9587920835627557e-08,
        "M_residual": 0.0,
        "minimum_psi": 0.17356235258058464,
        "relative_refinement": 9.92600152585192e-11
      },
      {
        "L": 5.5,
        "eta": 0.1,
        "H_residual": 2.5061230007894508e-08,
        "M_residual": 1.4898672406336324e-17,
        "minimum_psi": 0.17356235268479628,
        "relative_refinement": 9.925989214777441e-11
      },
      {
        "L": 5.5,
        "eta": 0.3,
        "H_residual": 2.912650615138901e-08,
        "M_residual": 4.4696017425804125e-17,
        "minimum_psi": 0.17356235351849642,
        "relative_refinement": 9.925755309726746e-11
      }
    ],
    "momentum_max": 3.043798496398225e-14,
    "H_offgrid_max": 5.445076944887717e-08,
    "boundary_max": 3.280948934570558e-12,
    "relative_refinement_max": 9.92600152585192e-11,
    "even_source_integral": -7.692045499692535e-06,
    "odd_source_integral": 0.0,
    "omitted_momentum_residual": 6.340254806093677e-07,
    "seam_max": 2.7755575615628914e-17,
    "missing_scalar_sign": 1.7317897831296805,
    "physical_H": [
      4.123368769178709e-06,
      3.9346830959896563e-07,
      5.480556209190615e-07
    ],
    "physical_M": [
      2.9023787225953904e-13,
      1.8112811105007222e-14,
      1.1442376679943052e-15
    ],
    "physical_ratios": [
      0.7179349952458096,
      15.829588215494203
    ],
    "wrong_f_H": 0.0714286353338755,
    "null_min": 2.4390094414377614e-05,
    "localization": [
      {
        "L": 3.5,
        "eta": 0.0,
        "bulk_error": 0.15040061226721269,
        "radii": [
          0.8495993877327873,
          0.2148056925779483,
          0.21444079355386245,
          0.2143192444643684
        ],
        "neck_ratio": 0.25225917951317456,
        "theta_plus": -1.3676207916587786e-10,
        "theta_minus": 1.3676207916587786e-10
      },
      {
        "L": 3.5,
        "eta": 0.1,
        "bulk_error": 0.1504006357436144,
        "radii": [
          0.8495993642563856,
          0.2148057105977928,
          0.21444081160536174,
          0.21431926252643302
        ],
        "neck_ratio": 0.25225920774318916,
        "theta_plus": 0.0006376153510711153,
        "theta_minus": 0.0006376152484223146
      },
      {
        "L": 3.5,
        "eta": 0.3,
        "bulk_error": 0.15040082355447892,
        "radii": [
          0.8495991764455211,
          0.21480585475624267,
          0.21444095601704952,
          0.21431940702261104
        ],
        "neck_ratio": 0.25225943358285946,
        "theta_plus": 0.001912841921221325,
        "theta_minus": 0.0019128421393110374
      },
      {
        "L": 4.5,
        "eta": 0.0,
        "bulk_error": 0.05359909685173048,
        "radii": [
          0.9464009031482695,
          0.08568333384027023,
          0.08552522532527378,
          0.08547256519482561
        ],
        "neck_ratio": 0.09031327517809322,
        "theta_plus": 7.419306709594768e-11,
        "theta_minus": -7.419306709594768e-11
      },
      {
        "L": 4.5,
        "eta": 0.1,
        "bulk_error": 0.05359909763012838,
        "radii": [
          0.9464009023698716,
          0.08568333464630111,
          0.08552522613209072,
          0.08547256600190407
        ],
        "neck_ratio": 0.09031327610516136,
        "theta_plus": 0.0005016552974527776,
        "theta_minus": 0.0005016553638920111
      },
      {
        "L": 4.5,
        "eta": 0.3,
        "bulk_error": 0.05359910385731004,
        "radii": [
          0.94640089614269,
          0.08568334109454669,
          0.0855252325866233,
          0.08547257245853715
        ],
        "neck_ratio": 0.09031328352171208,
        "theta_plus": 0.0015049656766952036,
        "theta_minus": 0.001504965625224736
      },
      {
        "L": 5.5,
        "eta": 0.0,
        "bulk_error": 0.019387472409855278,
        "radii": [
          0.9806125275901447,
          0.03228420553194528,
          0.032223883816700785,
          0.03220379325994891
        ],
        "neck_ratio": 0.03284048730143162,
        "theta_plus": 1.2850544668274102e-10,
        "theta_minus": -1.2850544668274102e-10
      },
      {
        "L": 5.5,
        "eta": 0.1,
        "bulk_error": 0.019387472444069576,
        "radii": [
          0.9806125275559304,
          0.03228420557057539,
          0.03222388385536257,
          0.03220379329862101
        ],
        "neck_ratio": 0.03284048734201413,
        "theta_plus": 0.0004671110094620103,
        "theta_minus": 0.0004671110875796634
      },
      {
        "L": 5.5,
        "eta": 0.3,
        "bulk_error": 0.01938747271778407,
        "radii": [
          0.9806125272822159,
          0.03228420587961595,
          0.03222388416465713,
          0.032203793608000394
        ],
        "neck_ratio": 0.03284048766667681,
        "theta_plus": 0.00140133306963863,
        "theta_minus": 0.0014013331407113888
      }
    ],
    "reversal_scalar_difference": 0.0,
    "reconstruction_knot_error": 3.2809489390842175e-12,
    "reconstruction_relative_change": 2.1687782569441334e-12
  },
  "verdicts": {
    "FOUR_SCALAR_HANDLE_CONSTRAINT_DATA": false,
    "LOCALIZED_BULK_MOUTH_INITIAL_DATA": false
  },
  "unestablished": [
    "traversability",
    "crossing_evolution",
    "momentum_transfer_events",
    "sector_selection",
    "discrete_action",
    "quantum_statistics"
  ],
  "original_gates": {
    "action": true,
    "momentum": true,
    "hamiltonian": false,
    "seam": true,
    "physical": false,
    "localization": true,
    "controls": true,
    "evidence": true
  },
  "refinement_freeze": "dbec68f2b6de8f745c269167a2ac7a41a38f4654"
}
```
