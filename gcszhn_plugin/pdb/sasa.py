from pymol import cmd

def set_sasa_color(
        selection:str = "(all)", 
        level: str = "A", 
        palette: str = 'red_white_blue',
        minimum: float = None,
        maximum: float = None):
    """
    Annotated color by sasa.

    Parameters
    ---------------
    - selection: str
        pymol selection str.
    - level: str
        'A' for atomic level,
        'R' for residual level,
        'C' for chain level.
    - palette: str
        color palette of pymol.
    
    """
    
    if type(minimum) is not type(maximum):
        raise ValueError("Please specific minimum and maximum both or not!")
    
    auto_bound = minimum is None
    
    if auto_bound:
        minimum = float('inf')
        maximum = -minimum

    dot_solvent = cmd.get('dot_solvent')
    sasa_buffer = dict()

    cmd.set('dot_solvent', 'on')
    cmd.get_area(selection, load_b=1)

    def set_sasa(resi, chain, b):
        if level == "R":
            res_sel_str = f'resi {resi} and chain {chain} and {selection}'
        elif level == "C":
            res_sel_str = f'chain {chain} and {selection}'

        if level == 'A':
            v = b
        else:
            if res_sel_str not in sasa_buffer:
                sasa_buffer[res_sel_str] = sum(atom.b for atom in cmd.get_model(res_sel_str).atom)
            v = sasa_buffer[res_sel_str]
        
        if auto_bound:
            nonlocal minimum, maximum
            minimum = min(minimum, v)
            maximum = max(maximum, v)
        return v
    
    cmd.alter(
        selection,
        f'b = set_sasa(resi, chain, b)',
        space={'set_sasa': set_sasa})
    cmd.spectrum(
        'b',
        palette=palette,
        selection=selection,
        minimum=minimum,
        maximum=maximum)

    cmd.set('dot_solvent', dot_solvent)