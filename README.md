# schrodinger_utils

Various integrations to the [Schrodinger](https://www.schrodinger.com/) ecosystem. With major interest on 
analysis tools for Molecular Dynamics simulations run with Desmond.

Please understand that this work is primarily meant to be used by co-workers and students of the LBBC lab @ University of Milan. Majority of the code was written during my Master/PhD and it was not originally meant to be shared, thus documentation must be improved.

Although -for the time being- this code is not optimized for general use, you can still find useful stuff here and most plugins/scripts should work.
For anyone interested, please contact me at: tommaso.laurenzi@unimi.it or open an issue.

**IMPORTANT**:
- Some plugins/scripts are still work in progress.
- Some plugins/scripts were written a long time ago, when I first started programming and I had limited knowledge of the Schrodinger API. Thus, despite best efforts, *may not* contain state-of-the-art code.
- `trjpeek.sh` and `mkhosts.sh` scripts contain hardcoded parameters that only make sense if you work with us.
- Some major tools, like `mi.py` for calculating Generalized correlations matrix, and `flow.py` to track the paths of certain molecules are also actively developed [here](https://github.com/uliano/xlence_scripts) and supporting other formats via MDAnalyis. These tools may soon find a new home.

## Overview
- Maestro plugins
- Docs
- Shell scripts
- Input templates

### Plugins for Maestro
- **maestro_applycscheme.py**: Map data to atom colors, requires python repl, exports `map2color()` function.
- **maestro_dummy_at_center_of_mass.py**: Add dummy atom at center of mass of current selection as a separate entry.
- **maestro_interacting.py**: Select interacting atoms between two groups or in current selection.
- **maestro_moi.py**: Draw principal moments of inertia of selected atoms.
- **maestro_notes.py**: Add and manage notes. Notes are added to *project table* entries as properties
- **maestro_pbitches.py**: Gain control over periodic box wrapping in the workspace. Make PBCs your ...
- **maestro_renumber.py**: Better residue renumbering GUI.
- **maestro_resdist.py**: Calculate per-residue distances between multiple versions of the same structure. Draws plots and consensus structure.
- **maestro_resdist_v2.py**: WIP, no GUI. As above, but structures need not be the same, read residue mapping from alignment.
- **maestro_select_by_b_factor.py**: Select residues with B factor above threshold within selection.
- **maestro_transform.py**: GUI to Rotate and Translate workspace structures. Transformations are applied relative to selected reference or to absolute coordinates.
- **map_interactions_to_structure.py**: Functions to parse `trj_interactions.py` output and map residue interactions to the workspace structure.
- **parse_prime_interaction_energies.py**: Functions to parse `prime_calcenergy` output and map interaction energies to workspace structure.

### Docs
- **schrodinger_tips.md**: Some general tips and obscure features.

### Scripts
- **mi.py**: Compute Mutal Information (Shannon Entropy / Generalized Correlations) matrix of selected atoms within Desmond Molecular Dynamics trajectory.
- **trj_flow.py**: Track trajectories of molecules that pass within a certain region during Molecular Dynamics.
- **trj_act.py**: Align, Center and Translate stuff in Desmond trajectories.
- **trj_boxinfo.py**: Report MD box info.
- **trj_conv.py**: Export Desmond trajectory and cms topology to pdb+xtc format.
- **trj_ene.py**: Report energies of and between selected atom groups as calculated by Desmond vrun.
- **trj_fcluster.py**: RMSD based clustering of Desmond trajectory.
- **trj_gcluster.py**: Implementation of gromos clustering method.
- **trj_helix_axis.py**: Report principal axis of α-helices.
- **trj_interactions.py**: Calculate interactions between selections in MD.
- **trj_log.py**: State and ETA of a running MD job.
- **trj_query.py**: Convert times and frame numbers.
- **trj_measure.py**: WIP. Various geometry analyses.
- **trj_periodic_shortest_distance.py**: Measure distances between selections across PBCs.
- **trj_rmsd.py**: Calculate RMSD and RMSF of selections.
- **trj_sasa.py**: Calculate SASA.
- **trj_ssp.py**: Calculate Secondary Structure information.
- **kerseq2hills.py**: Convert Desmond metadynamics kerseq file into plumed readable HILLS file.
- **biasf.py**: Convert between biasf and kT.
- **enhsamp_patched.py**: Patched version of enhsamp.py that returns a string compatible with msj/cfg files.
- **trjpeek.sh**: Get a copy of a running trajectory from its scratch directory
- **mkhosts.sh**: Bash utility to make user and version -specific schrodinger.hosts based on a template
- **XL_dock.py**: Toy tool to perform protein-protein docking of structures based on cross-link data.
- **schrun**: wrapper shell script to allow completions.
- **\_schrun**: WIP zsh completion function. Its generation could be automated by parsing help messages.

### Templates and obscure input files commands
- **minimize.msj**
- **simulate.msj**
- **system-builder_template.csb**
- **system-setup.msj**
- **schrodinger.hosts.template**

