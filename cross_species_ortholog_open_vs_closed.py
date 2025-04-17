import subprocess

def cross_species_ortholog_open_vs_closed (
        tissue: str,
        direction: str,
        species_from: str,
        species_to: str
    ):
    """
    Identify the open chromatin regions in a given species whose orthologs in the other species are open or closed.
    
    Args:
        tissue (str): 'liver' or 'pancreas'
        direction (str): 'HumanToMouse' or 'MouseToHuman'
        species_from (str): 'human' or 'mouse' (source species whose open chromatin peaks were lifted to the other species using halLifter and HALPER)
        species_to (str): 'mouse' or 'human' (target species used to determine whether the orthologous region is open or closed)
    """
    # Construct filenames
    halper_file = f"idr.optimal_peak_{tissue}.{direction}.HALPER.narrowPeak"
    cleaned_halper_file = halper_file + ".cleaned"
    native_file = f"{species_to}_{tissue}_idr.optimal_peak.narrowPeak.cleaned"

    # Step 1: Clean HALPER file (cut first 3 columns)
    print(f"\n[STEP 1] Cleaning {halper_file}")
    with open(cleaned_halper_file, "w") as out:
        subprocess.run(["cut", "-f1-3", halper_file], stdout=out)

    # Step 2: Intersect for conserved (ortholog open)
    conserved_file = f"{species_from}_to_{species_to}_{tissue}_conserved.bed"
    print(f"[STEP 2] Finding conserved peaks: {conserved_file}")
    with open(conserved_file, "w") as out:
        subprocess.run(["bedtools", "intersect", "-a", cleaned_halper_file, "-b", native_file, "-u"], stdout=out)

    # Step 3: Intersect for non-conserved (ortholog closed)
    closed_file = f"{species_from}_to_{species_to}_{tissue}_closed.bed"
    print(f"[STEP 3] Finding non-conserved peaks: {closed_file}")
    with open(closed_file, "w") as out:
        subprocess.run(["bedtools", "intersect", "-a", cleaned_halper_file, "-b", native_file, "-v"], stdout=out)

    # Step 4: Count peaks
    print("[STEP 4] Peak counts:")
    def count_lines(file): return int(subprocess.check_output(["wc", "-l", file]).decode().split()[0])
    
    lifted_total = count_lines(cleaned_halper_file)
    conserved_total = count_lines(conserved_file)
    closed_total = count_lines(closed_file)

    print(f"  Total lifted peaks:        {lifted_total}")
    print(f"  Conserved (open) peaks:    {conserved_total}")
    print(f"  Not conserved (closed):    {closed_total}")
