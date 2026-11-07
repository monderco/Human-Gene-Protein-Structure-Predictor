from pathlib import Path
import webbrowser

import py3Dmol
import requests


OUTPUT = Path("results")
OUTPUT.mkdir(exist_ok=True)


def get_json(url, params=None):
    response = requests.get(url, params=params, timeout=30)
    response.raise_for_status()
    return response.json()


def find_human_proteins(gene):
    """Find reviewed human proteins with this exact gene name."""
    data = get_json(
        "https://rest.uniprot.org/uniprotkb/search",
        params={
            "query": (
                f"(gene_exact:{gene}) AND "
                "(organism_id:9606) AND (reviewed:true)"
            ),
            "format": "json",
            "size": 25,
        },
    )

    ## Confirm that the submitted name matches a primary name or synonym ##
    matches = []

    for entry in data.get("results", []):
        names = []

        for gene_info in entry.get("genes", []):
            primary = gene_info.get("geneName", {}).get("value")
            if primary:
                names.append(primary)

            names.extend(
                synonym["value"]
                for synonym in gene_info.get("synonyms", [])
            )

        if gene.casefold() in {name.casefold() for name in names}:
            matches.append(entry)

    return matches


def choose_protein(entries):
    if not entries:
        raise ValueError(
            "No reviewed human protein found. Check the gene symbol."
        )

    if len(entries) == 1:
        return entries[0]

    print("\nMultiple proteins found:")

    for number, entry in enumerate(entries, start=1):
        description = entry.get("proteinDescription", {})
        name = (
            description.get("recommendedName", {})
            .get("fullName", {})
            .get("value", "Unnamed protein")
        )
        print(f"{number}. {entry['primaryAccession']} — {name}")

    choice = int(input("Choose a protein number: "))

    if not 1 <= choice <= len(entries):
        raise ValueError("Invalid protein number.")

    return entries[choice - 1]


def main():
    gene = input("Enter a human gene symbol: ").strip()

    # Keep the input suitable for a UniProt gene-symbol query.
    if not gene or not all(
        character.isalnum() or character in "-_."
        for character in gene
    ):
        raise ValueError("Enter a gene symbol using letters, numbers, -, _, or .")

    protein = choose_protein(find_human_proteins(gene))
    accession = protein["primaryAccession"]
    sequence = protein["sequence"]["value"]

    print(f"Found UniProt protein: {accession}")

    response = requests.get(
        f"https://alphafold.ebi.ac.uk/api/prediction/{accession}",
        timeout=30,
    )

    if response.status_code == 404:
        raise ValueError("No AlphaFold prediction is available for this protein.")

    response.raise_for_status()
    predictions = response.json()

    if not predictions:
        raise ValueError("No AlphaFold prediction is available for this protein.")

    ## This beginner version supports predictions covering the full sequence ##
    matches = [
        prediction
        for prediction in predictions
        if prediction.get("uniprotSequence") == sequence
        and prediction.get("pdbUrl")
    ]

    if not matches:
        raise ValueError(
            "No matching full-length prediction found. "
            "Partial structures need additional handling."
        )

    prediction = matches[0]

    response = requests.get(prediction["pdbUrl"], timeout=60)
    response.raise_for_status()
    pdb_text = response.text

    pdb_path = OUTPUT / f"{accession}.pdb"
    pdb_path.write_text(pdb_text, encoding="utf-8")

    fasta_path = OUTPUT / f"{accession}.fasta"
    fasta_path.write_text(
        f">{accession} gene={gene}\n{sequence}\n",
        encoding="utf-8",
    )

    viewer = py3Dmol.view(width=900, height=650)
    viewer.addModel(pdb_text, "pdb")
    viewer.setStyle({}, {"cartoon": {"color": "spectrum"}})
    viewer.setBackgroundColor("white")
    viewer.zoomTo()

    html_path = OUTPUT / f"{accession}.html"
    viewer.write_html(str(html_path))

    print(f"Saved protein sequence: {fasta_path}")
    print(f"Saved structure: {pdb_path}")
    print(f"Opening interactive structure: {html_path}")

    webbrowser.open(html_path.resolve().as_uri())


if __name__ == "__main__":
    try:
        main()
    except (requests.RequestException, ValueError) as error:
        print(f"Could not complete the request: {error}")