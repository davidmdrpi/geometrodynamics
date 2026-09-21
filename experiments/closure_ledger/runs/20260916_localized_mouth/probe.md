# Localized four-scalar handle experiment

Freeze: `66bce68ce5c9e07f4369283b18cb1ac4d3e94f73`.

Frozen gates: 6/8.

```json
{
  "gates": {
    "action": true,
    "momentum": true,
    "hamiltonian": false,
    "seam": true,
    "physical": false,
    "localization": true,
    "controls": true,
    "evidence": true
  },
  "passed": 6,
  "total": 8,
  "metrics": {
    "cases": [
      {
        "L": 3.5,
        "eta": 0.0,
        "H_residual": 4.856448044820993e-07,
        "M_residual": 0.0,
        "minimum_psi": 0.44774691070344674,
        "relative_refinement": 4.1178519677127935e-12
      },
      {
        "L": 3.5,
        "eta": 0.1,
        "H_residual": 4.855871809095191e-07,
        "M_residual": 1.0145994423305452e-14,
        "minimum_psi": 0.44774692957070217,
        "relative_refinement": 4.117696111337193e-12
      },
      {
        "L": 3.5,
        "eta": 0.3,
        "H_residual": 4.856383956641785e-07,
        "M_residual": 3.043798496398225e-14,
        "minimum_psi": 0.4477470805083964,
        "relative_refinement": 4.1176961968455176e-12
      },
      {
        "L": 4.5,
        "eta": 0.0,
        "H_residual": 3.8306762142892303e-07,
        "M_residual": 0.0,
        "minimum_psi": 0.2827583132852744,
        "relative_refinement": 8.787986359405928e-11
      },
      {
        "L": 4.5,
        "eta": 0.1,
        "H_residual": 3.8305837039842583e-07,
        "M_residual": 4.200630276569658e-16,
        "minimum_psi": 0.28275831462025475,
        "relative_refinement": 8.787986353831038e-11
      },
      {
        "L": 4.5,
        "eta": 0.3,
        "H_residual": 3.8298587791418814e-07,
        "M_residual": 1.2601890697360075e-15,
        "minimum_psi": 0.2827583253000938,
        "relative_refinement": 8.787965462389775e-11
      },
      {
        "L": 5.5,
        "eta": 0.0,
        "H_residual": 4.281440644238188e-07,
        "M_residual": 0.0,
        "minimum_psi": 0.17356235258058383,
        "relative_refinement": 9.925176716535392e-11
      },
      {
        "L": 5.5,
        "eta": 0.1,
        "H_residual": 4.2825092910758755e-07,
        "M_residual": 1.4898672406336324e-17,
        "minimum_psi": 0.17356235268479636,
        "relative_refinement": 9.925164405460953e-11
      },
      {
        "L": 5.5,
        "eta": 0.3,
        "H_residual": 4.2827458607286317e-07,
        "M_residual": 4.4696017425804125e-17,
        "minimum_psi": 0.17356235351849655,
        "relative_refinement": 9.92517671214675e-11
      }
    ],
    "momentum_max": 3.043798496398225e-14,
    "H_offgrid_max": 4.856448044820993e-07,
    "boundary_max": 3.2526065174565133e-19,
    "relative_refinement_max": 9.925176716535392e-11,
    "even_source_integral": -7.692045499692535e-06,
    "odd_source_integral": 0.0,
    "omitted_momentum_residual": 6.340254806093677e-07,
    "seam_max": 2.7755575615628914e-17,
    "missing_scalar_sign": 1.7317897831296805,
    "physical_H": [
      0.0010971068796524776,
      0.00027598283477686615,
      0.000171257720012433
    ],
    "physical_M": [
      1.5950848163191876e-10,
      3.9877104195773176e-11,
      9.969855141652341e-12
    ],
    "physical_ratios": [
      1.61150594996144,
      3.999767662538394
    ],
    "wrong_f_H": 0.07153088002237476,
    "null_min": 2.4390094414677562e-05,
    "localization": [
      {
        "L": 3.5,
        "eta": 0.0,
        "bulk_error": 0.15040061226721269,
        "radii": [
          0.8495993877327873,
          0.2148056925779481,
          0.2144407935538624,
          0.21431924446437028
        ],
        "neck_ratio": 0.2522591795131768,
        "theta_plus": 0.0,
        "theta_minus": 0.0
      },
      {
        "L": 3.5,
        "eta": 0.1,
        "bulk_error": 0.1504006357436144,
        "radii": [
          0.8495993642563856,
          0.21480571059779266,
          0.2144408116053617,
          0.21431926252643133
        ],
        "neck_ratio": 0.25225920774318716,
        "theta_plus": 0.0006376152997467312,
        "theta_minus": 0.0006376152997467291
      },
      {
        "L": 3.5,
        "eta": 0.3,
        "bulk_error": 0.15040082355447892,
        "radii": [
          0.8495991764455211,
          0.2148058547562425,
          0.21444095601704943,
          0.21431940702261265
        ],
        "neck_ratio": 0.2522594335828614,
        "theta_plus": 0.001912842030266136,
        "theta_minus": 0.0019128420302661382
      },
      {
        "L": 4.5,
        "eta": 0.0,
        "bulk_error": 0.05359909685173048,
        "radii": [
          0.9464009031482695,
          0.08568333384013485,
          0.08552522532516128,
          0.08547256519482481
        ],
        "neck_ratio": 0.09031327517809237,
        "theta_plus": -1.7944367045439945e-17,
        "theta_minus": 1.7944367045439945e-17
      },
      {
        "L": 4.5,
        "eta": 0.1,
        "bulk_error": 0.05359909763012838,
        "radii": [
          0.9464009023698716,
          0.0856833346461657,
          0.08552522613197806,
          0.08547256600190414
        ],
        "neck_ratio": 0.09031327610516143,
        "theta_plus": 0.0005016553306723932,
        "theta_minus": 0.0005016553306723932
      },
      {
        "L": 4.5,
        "eta": 0.3,
        "bulk_error": 0.05359910385731004,
        "radii": [
          0.94640089614269,
          0.08568334109441113,
          0.08552523258651074,
          0.08547257245853668
        ],
        "neck_ratio": 0.09031328352171158,
        "theta_plus": 0.0015049656509599945,
        "theta_minus": 0.0015049656509599945
      },
      {
        "L": 5.5,
        "eta": 0.0,
        "bulk_error": 0.019387472409855278,
        "radii": [
          0.9806125275901447,
          0.03228420553184434,
          0.03222388381658778,
          0.03220379325994861
        ],
        "neck_ratio": 0.03284048730143131,
        "theta_plus": 3.87951745928469e-17,
        "theta_minus": -3.87951745928469e-17
      },
      {
        "L": 5.5,
        "eta": 0.1,
        "bulk_error": 0.019387472444069576,
        "radii": [
          0.9806125275559304,
          0.032284205570474414,
          0.0322238838552496,
          0.03220379329862104
        ],
        "neck_ratio": 0.032840487342014155,
        "theta_plus": 0.00046711104852083555,
        "theta_minus": 0.00046711104852083555
      },
      {
        "L": 5.5,
        "eta": 0.3,
        "bulk_error": 0.01938747271778407,
        "radii": [
          0.9806125272822159,
          0.032284205879515006,
          0.032223884164544137,
          0.03220379360800044
        ],
        "neck_ratio": 0.03284048766667687,
        "theta_plus": 0.0014013331051750027,
        "theta_minus": 0.0014013331051750027
      }
    ],
    "reversal_scalar_difference": 0.0
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
  ]
}
```
