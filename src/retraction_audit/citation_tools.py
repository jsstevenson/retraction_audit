"""Tools for handling different kinds of citations.

In a lot of our use cases, it's easiest to work with PMIDs, so these methods help to
lift over other kinds of citations.
"""
import logging
import os
import sqlite3
from pathlib import Path
from typing import Optional

from Bio import Entrez

_logger = logging.getLogger(__name__)


Entrez.email = os.environ.get("NCBI_AUTH_EMAIL")


class SqlitePmidCache:
    """Cache PMID lookup results.

    Generally, DOI/PMCID to PMID references should be static and unchanging. Therefore,
    particularly during development/testing, we should reuse the values of those lookups
    between runs.
    """

    def __init__(self, db_file: Path) -> None:
        """Initialize cache instance

        :param db_file: path to sqlite DB file (doesn't need to exist yet)
        """
        self.con = sqlite3.connect(db_file)
        with self.con:
            cur = self.con.cursor()
            cur.execute(
                """
            CREATE TABLE IF NOT EXISTS doi_to_pmid(
                doi TEXT PRIMARY KEY,
                pmid TEXT
            );
            """
            )
            cur.execute(
                """
            CREATE TABLE IF NOT EXISTS pmcid_to_pmid(
                pmcid TEXT PRIMARY KEY,
                pmid TEXT
            );
            """
            )
            cur.execute(
                """
            CREATE TABLE IF NOT EXISTS author_title_to_pmid(
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                author TEXT,
                title TEXT,
                pmid TEXT
            )
            """
            )

    def _get_pmid_from_db(self, query: str, *args) -> str:
        cur = self.con.cursor()
        res = cur.execute(query, args)
        value = res.fetchone()
        if value is None:
            raise KeyError
        return value[0]

    def get_pmid_from_doi(self, doi: str) -> str:
        query = "SELECT * FROM doi_to_pmid WHERE doi = ?;"
        return self._get_pmid_from_db(query, doi)

    def get_pmid_from_pmcid(self, pmcid: str) -> str:
        query = "SELECT * FROM pmcid_to_pmid WHERE pmcid = ?;"
        return self._get_pmid_from_db(query, pmcid)

    def get_pmid_from_author_title(self, author: str, title: str) -> str:
        query = """SELECT * FROM author_title_to_pmid
        WHERE author COLLATE NOCASE = ? AND title COLLATE NOCASE = ?;
        """
        return self._get_pmid_from_db(query, author, title)

    def _add_pmid_to_db(self, query: str, *args) -> None:
        with self.con:
            cur = self.con.cursor()
            cur.execute(query, args)

    def add_pmid_from_doi(self, doi: str, pmid: str) -> None:
        self._add_pmid_to_db(
            "INSERT INTO doi_to_pmid(doi, pmid) VALUES(?, ?);",
            doi,
            pmid,
        )

    def add_pmid_from_pmcid(self, pmcid: str, pmid: str) -> None:
        self._add_pmid_to_db(
            "INSERT INTO pmcid_to_pmid(pmcid, pmid) VALUES(?, ?);",
            pmcid,
            pmid,
        )

    def add_pmid_from_author_title(self, author: str, title: str, pmid: str) -> None:
        self._add_pmid_to_db(
            "INSERT INTO author_title_to_pmid(author, title, pmid) VALUES (?, ?, ?);",
            author,
            title,
            pmid,
        )


def get_pmid_from_doi(doi: str, pmid_cache: Optional[SqlitePmidCache]) -> str:
    """Get PMID of article via DOI

    :param doi: article DOI
    :param pmid_cache:
    :return: PMID (as string)
    :raise KeyError: if no PMID available
    """
    if pmid_cache:
        try:
            return pmid_cache.get_pmid_from_doi(doi)
        except KeyError:
            pass
    with Entrez.esearch(db="pubmed", term=f"{doi}[DOI]") as h:
        result = Entrez.read(h)
    try:
        result = result["IdList"][0]
    except IndexError as e:
        raise KeyError from e
    if pmid_cache:
        pmid_cache.add_pmid_from_doi(doi, result)
    return result


def get_pmid_from_pmcid(pmcid: str, pmid_cache: Optional[SqlitePmidCache]) -> str:
    """Get PMID of article via PMC ID

    :param pmcid: article PMCID. Will remove leading ``"PMC"`` if included.
    :param pmid_cache:
    :return: PMID (as string)
    :raise KeyError: if no PMID available
    """
    pmcid = pmcid.upper()
    if pmcid.startswith("PMC"):
        pmcid = pmcid[3:]
    if pmid_cache:
        try:
            return pmid_cache.get_pmid_from_pmcid(pmcid)
        except KeyError:
            pass
    with Entrez.elink(dbfrom="pmc", id=pmcid) as h:
        result = Entrez.read(h)
    links = result[0]["LinkSetDb"]
    for link in links:
        if link.get("LinkName") == "pmc_pubmed":
            result = link["Link"][0]["Id"]
            break
    else:
        raise KeyError
    if pmid_cache:
        pmid_cache.add_pmid_from_pmcid(pmcid, result)
    return result


def get_pmid_from_author_title(
    author: str, title: str, pmid_cache: Optional[SqlitePmidCache]
) -> str:
    if pmid_cache:
        try:
            return pmid_cache.get_pmid_from_author_title(author, title)
        except KeyError:
            pass
    search_term = f"({author}[AUTHOR]) AND ({title}[TITLE])"
    with Entrez.esearch(db="pubmed", term=search_term) as h:
        result = Entrez.read(h)
    num_results = len(result.get("IdList", []))
    if num_results > 1:
        _logger.error(
            "Got unexpected number of PMID matches for %s by %s: %s",
            title,
            author,
            num_results,
        )
        raise KeyError
    if num_results == 0:
        _logger.debug("Got 0 PMID matches for %s by %s: %s", title, author, num_results)
        raise KeyError
    pmid = result["IdList"][0]
    if pmid_cache:
        pmid_cache.add_pmid_from_author_title(author, title, pmid)
    return pmid
