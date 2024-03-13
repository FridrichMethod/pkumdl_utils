import os
import gzip
import multiprocessing as mp
from typing import Iterator
from warnings import warn

import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import FilterCatalog
from rdkit.ML.Cluster import Butina
from rdkit.ML.Descriptors import MoleculeDescriptors
from PIL import Image

SIM_FUNCS = {name.lower(): func for name, func, _ in DataStructs.similarityFunctions}


def read_mols(file: str) -> Iterator[Chem.rdchem.Mol]:
    """Read molecules from a file

    Parameters
    ----------
    file : str
        Path to the file, which can be a .sdf

    Returns
    -------
    mol_gen : Iterator[rdkit.Chem.rdchem.Mol]
        A generator of molecules

    Notes
    -----
    Using multithreading and multiprocessing to read and process molecules from a file, so the order of the molecules may be different from the original file.
    """

    # Get the file extension
    ext = os.path.splitext(file)[1].lower()
    thread_num = 0
    queue_size = 1000
    match ext:
        case ".sdf":
            return (
                mol
                for mol in Chem.MultithreadedSDMolSupplier(
                    file,
                    numWriterThreads=thread_num,
                    sizeInputQueue=queue_size,
                    sizeOutputQueue=queue_size,
                )
                if mol is not None
            )
        case ".sdfgz":
            return (
                mol
                for mol in Chem.MultithreadedSDMolSupplier(
                    gzip.open(file),
                    numWriterThreads=thread_num,
                    sizeInputQueue=queue_size,
                    sizeOutputQueue=queue_size,
                )
                if mol is not None
            )
        case ".mae":
            return (mol for mol in Chem.MaeMolSupplier(file) if mol is not None)
        case ".maegz":
            return (
                mol for mol in Chem.MaeMolSupplier(gzip.open(file)) if mol is not None
            )
        case ".smi":
            return (
                mol
                for mol in Chem.MultithreadedSmilesMolSupplier(
                    file,
                    numWriterThreads=thread_num,
                    sizeInputQueue=queue_size,
                    sizeOutputQueue=queue_size,
                )
                if mol is not None
            )
        case ".csv" | ".xlsx" | ".xls":
            if ext == ".csv":
                df = pd.read_csv(file)
            else:
                df = pd.read_excel(file)

            def _process_df_row(row: pd.Series) -> Chem.rdchem.Mol:
                """Process a row of a dataframe to a molecule"""

                mol = Chem.MolFromSmiles(row["smiles"])
                mol.SetProp("_Name", row["title"])

                for prop, value in row.items():
                    if prop not in {"smiles", "title"}:
                        mol.SetProp(prop, str(value))

                return mol

            return (
                mol
                for _, row in df.iterrows()
                if (mol := _process_df_row(row)) is not None
            )
        case _:
            raise TypeError(
                "Should be a .sdf, .sdfgz, .mae, .maegz, .smi, .csv, .xlsx or .xls file."
            )


def is_pains(mol: Chem.rdchem.Mol) -> bool:
    """Identify and filiter out PAINS compounds

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        A molecule

    Returns
    -------
    pains : bool
        True if the molecule is a PAINS compound, False otherwise.
    """

    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog.FilterCatalog(params)
    return catalog.HasMatch(mol)


def calc_descs(
    mol: Chem.rdchem.Mol,
    desc_names: tuple[str, ...] = (
        "MolWt",
        "MolLogP",
        "NumHAcceptors",
        "NumHDonors",
        "FractionCSP3",
        "NumRotatableBonds",
        "RingCount",
        "TPSA",
        "qed",
    ),
) -> tuple[float, ...]:
    """Calculate molecular descriptors

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        A molecule
    desc_names : tuple, optional
        A list of descriptor names, by default (
            'MolWt', 'MolLogP', 'NumHAcceptors', 'NumHDonors',
            'FractionCSP3', 'NumRotatableBonds', 'RingCount', 'TPSA', 'qed',
        ).

    Returns
    -------
    descriptors : list
        A list of molecular descriptors.
        'MolWt', 'MolLogP', 'NumHAcceptors' and 'NumHDonors' are lipinski's rule of five descriptors.
        'FractionCSP3', 'NumRotatableBonds', 'RingCount', 'TPSA' and 'qed' are other common descriptors.
    """

    calc = MoleculeDescriptors.MolecularDescriptorCalculator(desc_names)
    return calc.CalcDescriptors(mol)


def draw_structures(
    mols: list[Chem.rdchem.Mol],
    *args: list[str],
    output_file: str | None = None,
    mols_per_row: int = 12,
    sub_img_size: tuple[float, float] = (300, 300),
) -> Image.Image | None:
    """Draw molecules

    Parameters
    ----------
    mols : list[rdkit.Chem.rdchem.Mol]
        A list of molecules
    *args : list[str]
        Lists of molecule properties, in the same order as mols
    output_file : str, optional
        Path to the output file, by default None
    mols_per_row : int, optional
        Number of molecules per row, by default 10
    sub_img_size : tuple[float, float], optional
        Size of each sub-image, by default (600, 600)

    Returns
    -------
    img : PIL.Image.Image | None
        A PIL image
    """

    Chem.rdDepictor.SetPreferCoordGen(True)

    legends = [" ".join(mol_prop) for mol_prop in zip(*args)]
    max_mols = len(mols)

    try:
        img = Draw.MolsToGridImage(
            mols,
            molsPerRow=mols_per_row,  # molsPerRow should not be too small
            subImgSize=sub_img_size,
            legends=legends,
            # returnPNG should be set EXPLICITLY to False
            # to avoid error in Jupyter Notebook
            returnPNG=False,
            maxMols=max_mols,
        )
    except RuntimeError:
        warn(
            "molsPerRow is too small, try to increase it.",
            RuntimeWarning,
        )
        return None
    else:
        if output_file:
            img.save(output_file)

        return img


def gen_fp(
    mol: Chem.rdchem.Mol,
    fp_type: str = "morgan",
) -> DataStructs.cDataStructs.ExplicitBitVect:
    """Generate molecular fingerprints

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        A molecule
    fp_type : str, optional
        Type of fingerprints, by default 'morgan'

    Returns
    -------
    fp : rdkit.DataStructs.cDataStructs.ExplicitBitVect
        Fingerprints
    """

    fp_type = fp_type.lower()

    match fp_type:
        case "morgan":
            return AllChem.GetMorganGenerator(radius=2, fpSize=1024).GetFingerprint(
                mol
            )  # ECFP4
        case "maccs":
            return Chem.rdMolDescriptors.GetMACCSKeysFingerprint(mol)
        case "rdkit":
            return AllChem.GetRDKitFPGenerator().GetFingerprint(mol)
        case "atom_pair":
            return AllChem.GetAtomPairGenerator().GetFingerprint(mol)
        case "topological_torsion":
            return AllChem.GetTopologicalTorsionGenerator().GetFingerprint(mol)
        case _:
            raise ValueError(
                "fp_type should be one of 'morgan', 'maccs', 'rdkit', 'atom_pair' and 'topological_torsion'."
            )


def calc_sim(
    fp1: DataStructs.cDataStructs.ExplicitBitVect,
    fp2: DataStructs.cDataStructs.ExplicitBitVect,
    similarity_metric: str = "tanimoto",
) -> float:
    """Calculate similarity between two fingerprints

    Parameters
    ----------
    fp1 : rdkit.DataStructs.cDataStructs.ExplicitBitVect
        A fingerprint
    fp2 : rdkit.DataStructs.cDataStructs.ExplicitBitVect
        A fingerprint
    similarity_metric : str, optional
        Similarity metric, by default 'tanimoto'

    Returns
    -------
    similarity : float
        Similarity between two fingerprints
    """

    if similarity_metric not in SIM_FUNCS:
        raise ValueError(
            "similarity_metric should be one of 'tanimoto', 'dice', 'cosine', 'sokal', 'russel', 'kulczynski', 'mcconnaughey', 'tversky' and 'asymmetric'."
        )

    return DataStructs.FingerprintSimilarity(
        fp1, fp2, metric=SIM_FUNCS[similarity_metric]
    )


def cluster_fps(
    fps: list[DataStructs.cDataStructs.ExplicitBitVect],
    cutoff: float = 0.6,
    similarity_metric: str = "tanimoto",
) -> list[tuple[int, ...]]:
    """Cluster the molecules by their Morgan fingerprints

    Parameters
    ----------
    fps : list[rdkit.DataStructs.cDataStructs.ExplicitBitVect]
        A list of Morgan fingerprints
    cutoff : float, optional
        Tanimoto similarity cutoff, by default 0.6

    Returns
    -------
    clusters : list[tuple[int, ...]]
        A list of clusters sorted by cluster size
    """

    nfps = len(fps)

    distances = [
        1 - calc_sim(fps[i], fps[j], similarity_metric=similarity_metric)
        for i in range(nfps)
        for j in range(i)
    ]  # Calculate the distance matrix

    clusters = list(
        Butina.ClusterData(distances, nfps, 1 - cutoff, isDistData=True)
    )  # Cluster the molecules by their Morgan fingerprints
    clusters.sort(
        key=len, reverse=True
    )  # Sort the clusters by the number of molecules in each cluster

    return clusters
