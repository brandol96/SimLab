# calculator definitions
def boolean_to_string(parameter):
    # convert logical variable into DFTB+ pattern
    if type(parameter == bool):
        if parameter:
            return 'Yes'
        else:
            return 'No'


def fetch_dftb_calc(cluster, **kwargs):
    from ase.calculators.dftb import Dftb
    label = kwargs.get('label')
    kpts = kwargs.get('kpts')
    max_force = kwargs.get('max_force')
    max_driver_steps = kwargs.get('max_driver_steps')
    lattice_opt = boolean_to_string(kwargs.get('lattice_opt'))
    fix_angles = boolean_to_string(kwargs.get('fix_angles'))
    fix_lengths = kwargs.get('fix_lengths')
    for i in range(2):
        fix_lengths[i]=boolean_to_string(fix_lengths[i])

    SCC = boolean_to_string(kwargs.get('SCC'))
    max_SCC = kwargs.get('max_SCC')
    max_SCC_steps = kwargs.get('max_SCC_steps')
    fermi_filling = kwargs.get('fermi_filling')
    use_LennardJones = kwargs.get('use_LennardJones')

    eVA_to_HaBohr = 0.01944689673

    if cluster:
        if use_LennardJones:
            calc = Dftb(label=label,
                        Driver_="ConjugateGradient",
                        Driver_MaxForceComponent=max_force * eVA_to_HaBohr,
                        Driver_MaxSteps=max_driver_steps,
                        Driver_MovedAtoms='1:-1',
                        Driver_AppendGeometries='Yes',
                        Hamiltonian_SCC=SCC,
                        Hamiltonian_SCCTolerance=max_SCC,
                        Hamiltonian_MaxSCCIterations=max_SCC_steps,
                        Hamiltonian_Filling=f"Fermi{{Temperature [K] = {fermi_filling} }}",
                        Hamiltonian_Dispersion='LennardJones{Parameters = UFFParameters{}}'
                        )
        else:
            calc = Dftb(label=label,
                        Driver_="ConjugateGradient",
                        Driver_MaxForceComponent=max_force * eVA_to_HaBohr,
                        Driver_MaxSteps=max_driver_steps,
                        Driver_MovedAtoms='1:-1',
                        Driver_AppendGeometries='Yes',
                        Hamiltonian_SCC=SCC,
                        Hamiltonian_SCCTolerance=max_SCC,
                        Hamiltonian_MaxSCCIterations=max_SCC_steps,
                        Hamiltonian_Filling=f"Fermi{{Temperature [K] = {fermi_filling} }}",
                        )
    else:
        if use_LennardJones:
            calc = Dftb(label=label,
                        kpts=kpts,
                        Driver_="ConjugateGradient",
                        Driver_MaxForceComponent=max_force * eVA_to_HaBohr,
                        Driver_MovedAtoms='1:-1',
                        Driver_LatticeOpt=lattice_opt,
                        Driver_FixAngles=fix_angles,
                        Driver_FixLengths=fix_lengths[0] + fix_lengths[1] + fix_lengths[2],
                        Driver_MaxSteps=max_driver_steps,
                        Driver_AppendGeometries='Yes',
                        Hamiltonian_SCC=SCC,
                        Hamiltonian_SCCTolerance=max_SCC,
                        Hamiltonian_MaxSCCIterations=max_SCC_steps,
                        Hamiltonian_Filling=f"Fermi{{Temperature [K] = {fermi_filling} }}",
                        Hamiltonian_Dispersion='LennardJones{Parameters = UFFParameters{}}'
                        )
        else:
            calc = Dftb(label=label,
                        kpts=kpts,
                        Driver_="ConjugateGradient",
                        Driver_MaxForceComponent=max_force * eVA_to_HaBohr,
                        Driver_MovedAtoms='1:-1',
                        Driver_LatticeOpt=lattice_opt,
                        Driver_FixAngles=fix_angles,
                        Driver_FixLengths=self.fixLengths[0] + self.fixLengths[1] + self.fixLengths[2],
                        Driver_MaxSteps=max_driver_steps,
                        Driver_AppendGeometries='Yes',
                        Hamiltonian_SCC=SCC,
                        Hamiltonian_SCCTolerance=max_SCC,
                        Hamiltonian_MaxSCCIterations=max_SCC_steps,
                        Hamiltonian_Filling=f"Fermi{{Temperature [K] = {fermi_filling} }}",
                        )
    return calc
