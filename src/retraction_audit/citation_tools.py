"""Tools for handling different kinds of citations.

In a lot of our use cases, it's easiest to work with PMIDs, so these methods help to
lift over other kinds of citations.
"""

import os
from functools import lru_cache

from Bio import Entrez


@lru_cache
def get_pmid_from_doi(doi: str) -> str:
    """Get PMID of article via DOI

    :param doi: article DOI
    :return: PMID (as string)
    :raise KeyError: if no PMID available
    """
    Entrez.email = os.environ.get("NCBI_AUTH_EMAIL")
    with Entrez.esearch(db="pubmed", term=f"{doi}[DOI]") as h:
        result = Entrez.read(h)
    try:
        return result["IdList"][0]
    except IndexError as e:
        raise KeyError from e


@lru_cache
def get_pmid_from_pmcid(pmcid: str) -> str:
    """Get PMID of article via PMC ID

    :param pmcid: article PMCID. Will add leading ``"PMC"`` if not provided.
    :return: PMID (as string)
    :raise KeyError: if no PMID available
    """
    pmcid = pmcid.upper()
    if pmcid.startswith("PMC"):
        pmcid = pmcid[3:]
    Entrez.email = os.environ.get("NCBI_AUTH_EMAIL")
    with Entrez.elink(dbfrom="pmc", id=pmcid) as h:
        result = Entrez.read(h)
    links = result[0]["LinkSetDb"]
    for link in links:
        if link.get("LinkName") == "pmc_pubmed":
            return link["Link"][0]["Id"]
    raise KeyError
