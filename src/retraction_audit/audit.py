"""Provide publication retraction information."""
import csv
from collections import namedtuple
from enum import StrEnum
from typing import ClassVar, Dict, List, Optional

from retraction_audit.fetch_data import get_latest_data


class FlagType(StrEnum):
    """Define types of retraction flags."""

    RETRACTED = "Retraction"
    CONCERNED = "Expression of concern"


RetractionRecord = namedtuple(
    "RetractionData",
    [
        "pmid",
        "doi",
        "article_title",
        "journal",
        "authors",
        "retraction_watch_refs",
        "retraction_doi",
        "retraction_pmid",
        "retraction_reasons",
        "retraction_type",
    ],
)


class Audit:
    """Handle retraction checks."""

    data: ClassVar[Dict] = {}

    def __init__(self) -> None:
        """Initialize Audit instance."""
        retraction_file = get_latest_data()

        with retraction_file.open() as f:
            reader = csv.DictReader(f)
            for row in reader:
                self._add_to_map(row)

    @staticmethod
    def _process_list_row(item: str, strip: str = "") -> List[str]:
        item_list = item.split(";")
        item_list = [i.strip(strip) for i in item_list]
        return list(filter(None, item_list))

    def _add_to_map(self, row: Dict) -> None:
        pmid = row["OriginalPaperPubMedID"]
        doi = row["OriginalPaperDOI"]
        if not pmid and not doi:
            return
        nature = row["RetractionNature"]
        if nature != FlagType.CONCERNED and nature != FlagType.RETRACTED:
            return
        retraction_record = RetractionRecord(
            pmid=pmid,
            doi=doi,
            article_title=row["Title"],
            journal=row["Journal"],
            authors=self._process_list_row(row["Author"]),
            retraction_watch_refs=self._process_list_row(row["URLS"]),
            retraction_doi=row["RetractionDOI"],
            retraction_pmid=row["RetractionPubMedID"],
            retraction_reasons=self._process_list_row(row["Reason"], "+"),
            retraction_type=FlagType.CONCERNED
            if nature == FlagType.CONCERNED
            else FlagType.RETRACTED,
        )
        if pmid and pmid != "Unavailable":
            key = pmid.lower()
        elif doi and doi != "Unavailable":
            key = doi.lower()
        else:
            return
        self.data[key] = retraction_record

    def get_retraction_by_id(self, identifier: str) -> Optional[RetractionRecord]:
        """Look up known retraction status given a PMID.

        :param identifier: article ID
        :return: retraction description, if available
        """
        return self.data.get(identifier.lower())
