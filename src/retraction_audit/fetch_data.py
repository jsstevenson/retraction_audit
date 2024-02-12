"""Acquire and manage retraction data."""
import datetime
import os
from pathlib import Path

from wags_tails import CustomData
from wags_tails.utils.downloads import download_http
from wags_tails.utils.versioning import DATE_VERSION_PATTERN


def get_latest_data(silent: bool = False) -> Path:
    """TODO"""

    def _latest_version_cb() -> str:
        """Get latest version value.

        The CrossRef-hosted Retraction Watch data is theoretically updated daily,
        but the specific time isn't posted, so we can't say for sure whether a new
        version has been posted on a given date. The simplest solution is to just
        check no more than once daily.

        :return: the latest feasible version (ie, the current date)
        """
        return datetime.datetime.now(tz=datetime.timezone.utc).strftime(
            DATE_VERSION_PATTERN
        )

    def _download_cb(version: str, outfile_path: Path) -> None:
        """TODO"""
        tqdm_params = {
            "disable": silent,
            "unit": "B",
            "ncols": 80,
            "unit_divisor": 1024,
            "unit_scale": True,
        }
        credentials = os.environ.get("RETRACTION_WATCH_EMAIL")
        if not credentials:
            msg = "Must provide email address under environment variable $RETRACTION_WATCH_EMAIL."
            raise ValueError(msg)
        download_http(
            f"https://api.labs.crossref.org/data/retractionwatch?{credentials}",
            outfile_path,
            tqdm_params=tqdm_params,
        )

    cd = CustomData(
        "retraction_watch",
        "csv",
        _latest_version_cb,
        _download_cb,
    )

    return cd.get_latest()[0]
