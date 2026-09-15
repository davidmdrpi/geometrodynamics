# Compact twisted-handle constraint experiment

Public freeze: `2b1bfca7b2650db7f8dcf74ca1c54ecd8bdffa17`.

Frozen gates: **8/8**.

| Gate | Pass |
|---|---|
| symbolic | True |
| momentum_completion | True |
| tensor_descent | True |
| hamiltonian | True |
| physical_constraints | True |
| minimal_section | True |
| controls | True |
| evidence | True |

```json
{
  "gates": {
    "symbolic": true,
    "momentum_completion": true,
    "tensor_descent": true,
    "hamiltonian": true,
    "physical_constraints": true,
    "minimal_section": true,
    "controls": true,
    "evidence": true
  },
  "passed": 8,
  "total": 8,
  "metrics": {
    "momentum_errors": [
      2.2209284968002785e-05,
      5.556888485289124e-06,
      1.389507680304991e-06
    ],
    "momentum_ratios": [
      3.9967123736238226,
      3.9991779563747434
    ],
    "ongrid_max": 3.009092974792793e-12,
    "offgrid_max": 3.127859082852069e-12,
    "grid_differences": [
      [
        9.992007221626409e-16,
        2.4424906541753444e-15
      ],
      [
        2.6645352591003757e-15,
        1.7763568394002505e-15
      ],
      [
        8.01581023779363e-14,
        2.886579864025407e-15
      ],
      [
        4.907407813448117e-12,
        1.9984014443252818e-15
      ],
      [
        2.8021429621105653e-10,
        2.6645352591003757e-15
      ]
    ],
    "minimum_psi": 1.0,
    "maximum_K_norm_on_frozen_grids": 1.23386314901813,
    "seam_max": 5.551115123125783e-12,
    "bad_seam_max": 0.3470483955074361,
    "physical_H": [
      2.838722140799236e-06,
      7.132438140158468e-07,
      1.8428338110560642e-07
    ],
    "physical_M": [
      8.225125726488322e-08,
      2.056282586801798e-08,
      5.140672411790227e-09
    ],
    "wrong_weight_H": 0.021294504367012836,
    "wrong_weight_M": 0.009531078554682642,
    "physical_sample_max_K_norm": 1.2334921110291315,
    "sections": [
      {
        "epsilon": 0.02,
        "areas": [
          12.567528244263842,
          12.56752963524661,
          12.567533804717918,
          12.56975422343707
        ],
        "expansions_min": 0.9998601905604787,
        "expansions_max": 0.9998657635167806
      },
      {
        "epsilon": 0.05,
        "areas": [
          12.573604441862944,
          12.573613130283603,
          12.57363917381856,
          12.587506593720168
        ],
        "expansions_min": 0.9991268885548055,
        "expansions_max": 0.9991616706604031
      },
      {
        "epsilon": 0.1,
        "areas": [
          12.595286726690519,
          12.595321406278252,
          12.595425358195987,
          12.650750589106199
        ],
        "expansions_min": 0.9965174602288541,
        "expansions_max": 0.99665589338363
      },
      {
        "epsilon": 0.2,
        "areas": [
          12.68174130333363,
          12.681878866335445,
          12.682291208977684,
          12.901354406444256
        ],
        "expansions_min": 0.9862247694664251,
        "expansions_max": 0.9867675557195436
      }
    ],
    "controls": {
      "amplitude_reversal": {
        "pde": 1.5649981310872363e-13,
        "psi_difference": 0.0
      },
      "time_reversal": {
        "pde": 1.5649981310872363e-13,
        "psi_difference": 0.0
      },
      "untwisted": {
        "pde": 1.6829593274536592e-13
      }
    }
  },
  "verdicts": {
    "COMPACT_TWISTED_VACUUM_CONSTRAINT_DATA": true,
    "MINIMAL_SECTION_WITH_NONTRIVIAL_MOMENTUM_COMPLETION": true
  },
  "unestablished": {
    "round_s3_embedding": true,
    "four_scalar_interface": true,
    "traversability": true,
    "mouth_evolution": true,
    "radiative_momentum_transfer": true,
    "crossing_events": true,
    "discrete_action": true,
    "quantum_statistics": true
  }
}
```

Initial constraints only. No traversability, crossing, momentum-transfer or quantum claim.
