"""Reproducible source provenance for a manually interpreted model audit.

Hashes and anchors make the reading checkable; they do not prove absence of
every possible interface. The result concerns the named baseline models.
"""
import ast
import hashlib
import subprocess
from pathlib import Path
from functools import lru_cache
import copy
from geometrodynamics.waves.field_apparatus import BASELINE

ROOT = Path(__file__).resolve().parents[2]
# path, symbols, dimension, boundary/geometry, source law, missing connection.
SPECS = (
 ('geometrodynamics/waves/nonlinear_supported_tt.py', ['conformal_rhs','constraints','initial_data'], '3+1',
  'Homogeneous compact S3; no excised boundary.',
  'A,q,M,L coupled; F=1/kappa-Q/(6 A^2).',
  'No mouth embedding, fiber restriction or surface action.'),
 ('geometrodynamics/waves/coupled_multiplet_response.py', ['CoupledSupport','full_geometry'], '3+1',
  'Homogeneous tensor perturbations of round supported ESU.',
  'Full improved stress varies with metric despite vanishing field perturbations.',
  'No localized mouth coordinate or interface law.'),
 ('geometrodynamics/waves/reciprocal_scalar_tt.py', ['quaternion_derivatives','HarmonicMultiplet'], '3+1 projection',
  'Five homogeneous tensors and scalar harmonics.',
  'Reciprocal scalar/TT coupling; not all Einstein constraints.',
  'No shape variable for the throat fiber.'),
 ('geometrodynamics/waves/physical_throat.py', ['VacuumThroat','profile_scalar_curvature'], '3-space initial data',
  'Scalar-flat spherical neck, time-symmetric slice and spatial C1 gluing.',
  "fprime^2=1-f0/f; ambient matching fixes f0=sin(a)^3.",
  'No lapse/evolution or four-scalar perturbative junction.'),
 ('geometrodynamics/waves/areal.py', ['TubeModel','mouth_sphere','solve_matching'], '3-space constraints',
  'S2 mouths and chosen round product tube; monopole/dipole boundary channels.',
  'Constraint response and Dirichlet-to-Neumann matching.',
  'Local S2 channels are not the throat winding fiber; no m=2 shape evolution.'),
 ('geometrodynamics/waves/throat_operator.py', ['BoundaryCondition','MouthPair','mouth_flux'], 'S3 point extensions',
  'Hermitian self-adjoint boundary family on punctured S3.',
  'Im(q* phi_regular) flux; boundary matrix is supplied.',
  'Conservation constrains the family but selects no four-field mouth law.'),
 ('geometrodynamics/waves/finite_throat.py', ['FiniteThroat','dtn_matrix','measure_the_enlarged_system_is_conservative'], 'ambient plus interval model',
  'Chosen tube area/length with reciprocal endpoint matching.',
  'Conservative enlarged ambient-plus-tube system; frequency-dependent elimination.',
  'Conserving wave boundary model, not an Einstein-scalar interface or shape equation.'),
 ('geometrodynamics/tangherlini/traversable_throat.py', ['throat_radius','stress_tensor','scattering_matrix'], '4+1',
  'Ultrastatic ds^2=-dt^2+ds^2+(s^2+a^2)dOmega3^2 benchmark.',
  'rho=0, negative radial pressure; supporting stress specified by the geometry.',
  'Required 5D stress not realized by the 4D quartet; no dimensional reduction.'),
 ('geometrodynamics/shells/junction.py', ['surface_stress','Region','Gluing'], 'D-dimensional Einstein shells',
  'Spherical thin shells and vacuum-region metric family.',
  'Israel jump uses constant Einstein coupling and prescribed surface equation of state.',
  'Neither nonspherical m=2 response nor matching for nonminimal F R action supplied.'),
 ('geometrodynamics/transaction/network.py', ['MouthPort','closure_offset','derived_loop_eigenvalue'], 'transport/network',
  'Supplied port times, clock offsets and transmission data.',
  'Clock/phase closure combines specified propagators.',
  'Does not determine local metric deformation or supporting stress.'),
 ('experiments/closure_ledger/throat_order_field_probe.py', ['landau_V','vortex_profile'], 'effective GL model',
  'Prescribed core and asymptotic order-parameter conditions.',
  'V=(lambda/4)(|q|^2-q0^2)^2 with independent coefficients.',
  'Order field added; microscopic potential and self-gravitating coupling underived.'),
 ('experiments/closure_ledger/throat_action_derivation_probe.py', ['S_energy','S_Hopf','antipodally_symmetric_partition'], 'reduced history ansatz',
  'Symmetric two-segment closed orbit.',
  'Closure action quantum supplied; partition imposed by symmetric ansatz.',
  'Not a local four-scalar surface action; cannot derive its boundary force.'),
 ('experiments/closure_ledger/minimal_mixing_interaction_probe.py', ['chi_harmonic','restrict_to_qubit','charge_conserving_extension'], 'x/chi lattice',
  'Chosen mouth envelope and cos(m chi+phase) potential.',
  'Delta k=2 potential plus separately supplied carrier interaction.',
  'Operator and carrier coupling not produced by supported Einstein evolution.'),
 ('experiments/closure_ledger/throat_apparatus_pointer_probe.py', ['mouth_winding_matrix','joint_pointer_probs'], 'quantum lattice apparatus',
  'Prepared carrier and complete winding pointer branches.',
  'joint_pointer_probs explicitly computes abs(state)^2.',
  'Quantum probabilities and prepared occupation do not derive classical event frequencies.'),
 ('experiments/closure_ledger/sigma_z_readout_capstone_probe.py', ['fiber_H','qubit_restriction'], 'Nchi=8 weighted graph',
  'R_mid prescribed; hopping 1/R_mid^2; flat kinetic/site measure.',
  'Ellipticity changes a material/graph operator.',
  'No covariant dimensional reduction establishes that hopping and kinetic measure.'),
 # Direct dependencies and the additional candidate returned by the search.
 ('geometrodynamics/waves/scalar_esu_support.py', ['improved_stress'], '3+1',
  'Pointwise smooth scalar geometry.', 'Improved conformal stress includes G_ab phi^2 and second derivatives.',
  'No nonspherical interface action or surface scalar conditions.'),
 ('geometrodynamics/waves/initial_data.py', ['constraint_operator_eigenvalue','regularised_green'], '3-space constraints',
  'Linear conformal initial-data response on round S3 with excisions.',
  'Elliptic operator has a degree-1 kernel.',
  'Elliptic means a PDE type here, not an elliptic-mouth mixer; no evolution.'),
 ('experiments/closure_ledger/configuration_space_emergence_probe.py', ['mode'], 'x/chi lattice',
  'Prescribed bridge links and winding modes.', 'Kinematic lattice embedding plus imported Bell-state machinery.',
  'No map from the quaternionic four-scalar Einstein solution to this lattice.'),
 ('geometrodynamics/shells/multipole.py', ['mutual_stiffness'], 'static Newtonian D-dimensional shells',
  'Concentric prescribed shell deformations; Laplace Green kernel.',
  'A genuine ell=2 mutual stiffness exists; shear modulus remains an input.',
  'No relativistic nonminimal-scalar interface, kinetic shape equation or fiber map.'),
 ('geometrodynamics/bulk/closure_grounding.py', ['tube_frame_energy'], 'conditional scalar tube',
  'Nine unconstrained matrix channels with fixed endpoints I and Ad_G.',
  'Gradient-energy minimization gives the desired restoring energy conditionally.',
  'Additional channels and endpoint matching; no identification with the quartet.'),
 ('geometrodynamics/transaction/derived_network.py', ['derived_throat'], '4+1 benchmark transport',
  'Imports traversable_throat scattering into network ports.',
  'Uses the supported benchmark with specified clock offset.',
  'Transport closure does not supply the benchmark stress or an m=2 shape source.'),
 ('geometrodynamics/waves/neck.py', ['NeckThroat','rayleigh_quotient'], 'excised 3-space plus tube',
  'Round balls removed; conservative field matching to a prescribed tube.',
  'Positive scalar quadratic form, higher multipoles decouple from one-channel tube.',
  'Fixed geometry and scalar wave action, not coupled Einstein-scalar mouth dynamics.'),
)

SEARCHES = (
 ('nonlinear_supported_tt|coupled_multiplet_response', ['geometrodynamics','experiments/closure_ledger/*.py']),
 ('elliptic|quadrupol|winding.?2', ['geometrodynamics/waves','geometrodynamics/shells','geometrodynamics/tangherlini','geometrodynamics/transaction']),
 ('surface_action|junction|matching', ['geometrodynamics','experiments/closure_ledger/*throat*probe.py']),
 ('Ginzburg|order.parameter', ['geometrodynamics','experiments/closure_ledger/*throat*probe.py']),
)


def git(*args):
    return subprocess.check_output(['git',*args], cwd=ROOT)


@lru_cache(None)
def _inventory():
    rows = []
    for path, symbols, dimension, boundary, law, gap in SPECS:
        raw = git('show', BASELINE+':'+path)
        tree = ast.parse(raw.decode())
        nodes = {n.name:n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
        anchors = []
        for symbol in symbols:
            node = nodes[symbol]
            snippet = '\n'.join(raw.decode().splitlines()[node.lineno-1:node.end_lineno])
            anchors.append(dict(symbol=symbol, line=node.lineno, sha256=hashlib.sha256(snippet.encode()).hexdigest()))
        rows.append(dict(path=path, commit=BASELINE, sha256=hashlib.sha256(raw).hexdigest(),
                         anchors=anchors, dimension=dimension, boundary=boundary, defining_law=law,
                         missing_connection=gap))
    searches = []
    for pattern, paths in SEARCHES:
        args = ['git','grep','-l','-I','-E',pattern,BASELINE,'--',*paths]
        run = subprocess.run(args,cwd=ROOT,capture_output=True,text=True)
        if run.returncode not in (0,1):
            raise RuntimeError(run.stderr)
        searches.append(dict(command=args, exit_code=run.returncode, output=run.stdout.splitlines()))
    return dict(baseline=BASELINE, rows=rows, searches=searches,
                interpretation='manual_source_audit_of_named_models_not_a_nonexistence_theorem',
                map=None, interface_action=None, physical_response=None,
                missing=['physical mouth and fiber embedding for the quartet',
                         'metric and scalar interface conditions for F R coupling',
                         'mouth shape action and reciprocal forcing',
                         'preparation and physical record law'])


def inventory():
    return copy.deepcopy(_inventory())


def valid(data):
    # An edited interpretation cannot certify itself with an intact status flag.
    return data == inventory()
