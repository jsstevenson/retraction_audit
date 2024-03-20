"""Provide tools for analysis notebooks."""

import datetime
import json
import os
import sqlite3
from typing import List, Tuple

from retraction_audit.audit import FlagType, RetractionLookup, RetractionRecord

os.environ["RETRACTION_WATCH_EMAIL"] = "james.sharpsteen@gmail.com"

lookup = RetractionLookup(True)

RetractionExperimentResult = List[RetractionRecord]


def _adapt_datetime(dt: datetime.datetime) -> str:
    return dt.isoformat()


def _convert_datetime(raw_dt: bytes) -> datetime.datetime:
    return datetime.datetime.fromisoformat(raw_dt.decode())


sqlite3.register_adapter(datetime.datetime, _adapt_datetime)
sqlite3.register_converter("datetime", _convert_datetime)


def _adapt_retraction_experiment_result(results: RetractionExperimentResult) -> str:
    transformed_results = []
    for result in results:
        data = result._asdict()
        data["retraction_type"] = str(data["retraction_type"])
        transformed_results.append(data)
    return json.dumps(transformed_results)


def _convert_retraction_experiment_result(
    raw_result: bytes,
) -> RetractionExperimentResult:
    deserialized = json.loads(raw_result.decode())
    results = []
    for raw_result in deserialized:
        raw_result["retraction_type"] = FlagType(raw_result["retraction_type"])
        record = RetractionRecord(**raw_result)
        results.append((raw_result[0], record))
    return results


sqlite3.register_converter("RETRACTION_BLOB", _convert_retraction_experiment_result)


def store_results(
    src_name: str,
    retraction_count: int,
    retraction_data: RetractionExperimentResult,
    total_citations: int,
) -> None:
    con = sqlite3.connect("experiment_results.db")
    with con:
        cur = con.cursor()

        table_query = """CREATE TABLE IF NOT EXISTS experiment_results (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            datetime DATETIME,
            src_name TEXT,
            retraction_count INTEGER,
            total_citations INTEGER,
            retraction_data RETRACTION_BLOB
        );
        """
        cur.execute(table_query)

        insert_query = """INSERT INTO experiment_results(
            datetime,
            src_name,
            retraction_count,
            total_citations,
            retraction_data
        )
        VALUES (?, ?, ?, ?, ?);
        """
        cur.execute(
            insert_query,
            (
                datetime.datetime.now(tz=datetime.timezone.utc),
                src_name,
                retraction_count,
                total_citations,
                _adapt_retraction_experiment_result(retraction_data),
            ),
        )


def get_all_results() -> List[Tuple]:
    con = sqlite3.connect("experiment_results.db")
    with con:
        cur = con.cursor()
        return cur.execute("SELECT * FROM experiment_results").fetchall()


def get_latest_results() -> List[Tuple]:
    con = sqlite3.connect("experiment_results.db")
    with con:
        cur = con.cursor()
        return cur.execute(
            """
            SELECT
                r.src_name,
                r.retraction_count,
                r.total_citations,
                r.datetime
            FROM experiment_results r
            JOIN (
                SELECT src_name, MAX(datetime) AS max_date
                FROM experiment_results
                GROUP BY src_name
            ) AS sub
            ON r.src_name = sub.src_name AND r.datetime = sub.max_date;
        """
        ).fetchall()
