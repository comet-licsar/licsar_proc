#!/usr/bin/env python3
"""Download Sentinel-1 AUX files needed by GAFA.

GAFA normally needs AUX_INS and AUX_CAL SAFE directories for the RAW
acquisition date. AUX_PP1 is only needed when raw.apply_b2b_gain is enabled.

Script by Muhammet Nergizci (RIP 2026)
"""

from __future__ import annotations

import argparse
import html
import re
import shutil
import sys
import tempfile
import zipfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin
from urllib.request import Request, urlopen


MPC_ROOT = "https://sar-mpc.eu"
AUX_PAGES = {
    "INS": "https://sar-mpc.eu/adf/aux_ins/",
    "CAL": "https://sar-mpc.eu/adf/aux_cal/",
    "PP1": "https://sar-mpc.eu/adf/aux_pp1/",
}


@dataclass(frozen=True)
class AuxProduct:
    name: str
    aux_type: str
    platform: str
    validity_start: datetime
    generation_time: datetime
    download_url: str


def parse_s1_time(value: str) -> datetime:
    return datetime.strptime(value, "%Y%m%dT%H%M%S")


def parse_product_or_date(
    product: str | None, platform: str | None, date: str | None
) -> tuple[str, datetime]:
    if product:
        name = Path(product).name
        match = re.search(r"^(S1[ABCD])_.*_(\d{8}T\d{6})_(\d{8}T\d{6})_", name)
        if not match:
            raise ValueError(f"Could not parse Sentinel-1 platform/date from product name: {name}")
        return match.group(1), parse_s1_time(match.group(2))

    if not date:
        raise ValueError("Provide either --product or both --platform and --date.")

    platform_match = re.fullmatch(r"S1[ABCD]", platform or "")
    if not platform_match:
        raise ValueError("--platform must be one of S1A, S1B, S1C, or S1D.")
    return platform, parse_s1_time(normalize_date(date))


def normalize_date(value: str) -> str:
    value = value.strip()
    if re.fullmatch(r"\d{8}", value):
        return f"{value}T000000"
    if re.fullmatch(r"\d{8}T\d{6}", value):
        return value
    raise ValueError("--date must be YYYYMMDD or YYYYMMDDTHHMMSS.")


def fetch_text(url: str) -> str:
    req = Request(url, headers={"User-Agent": "gafa-s1-aux-downloader/1.0"})
    with urlopen(req, timeout=60) as response:
        return response.read().decode("utf-8", errors="replace")


def discover_last_page(page_text: str) -> int:
    pages = [int(p) for p in re.findall(r"[?&]page=(\d+)", page_text)]
    return max(pages, default=1)


def parse_aux_products(page_text: str, aux_type: str) -> list[AuxProduct]:
    products: list[AuxProduct] = []
    pattern = re.compile(
        r'href="(?P<href>/download/[^"]+)">\s*(?P<name>S1[ABCD]_AUX_'
        + re.escape(aux_type)
        + r"_V(?P<valid>\d{8}T\d{6})_G(?P<generated>\d{8}T\d{6}))",
        re.IGNORECASE,
    )
    for match in pattern.finditer(page_text):
        name = html.unescape(match.group("name"))
        products.append(
            AuxProduct(
                name=name,
                aux_type=aux_type,
                platform=name[:3],
                validity_start=parse_s1_time(match.group("valid")),
                generation_time=parse_s1_time(match.group("generated")),
                download_url=urljoin(MPC_ROOT, html.unescape(match.group("href"))),
            )
        )
    return products


def list_aux_products(aux_type: str) -> list[AuxProduct]:
    base_url = AUX_PAGES[aux_type]
    first_page = fetch_text(base_url)
    products = parse_aux_products(first_page, aux_type)
    for page in range(2, discover_last_page(first_page) + 1):
        products.extend(parse_aux_products(fetch_text(f"{base_url}?page={page}"), aux_type))
    return products


def choose_aux_product(
    products: list[AuxProduct], platform: str, acquisition_time: datetime, aux_type: str
) -> AuxProduct:
    valid = [
        product
        for product in products
        if product.platform == platform and product.validity_start <= acquisition_time
    ]
    if not valid:
        raise FileNotFoundError(
            f"No {platform} AUX_{aux_type} product with validity_start <= "
            f"{acquisition_time:%Y%m%dT%H%M%S}."
        )
    return max(valid, key=lambda product: (product.validity_start, product.generation_time))


def download_file(url: str, out_file: Path, overwrite: bool = False) -> None:
    if out_file.exists() and out_file.stat().st_size > 1024 and not overwrite:
        print(f"Already exists: {out_file}")
        return

    out_file.parent.mkdir(parents=True, exist_ok=True)
    tmp_file = out_file.with_suffix(out_file.suffix + ".tmp")
    req = Request(url, headers={"User-Agent": "gafa-s1-aux-downloader/1.0"})
    try:
        with urlopen(req, timeout=120) as response, tmp_file.open("wb") as dst:
            shutil.copyfileobj(response, dst, length=1024 * 1024)
    except (HTTPError, URLError) as exc:
        if tmp_file.exists():
            tmp_file.unlink()
        raise RuntimeError(f"Download failed for {url}: {exc}") from exc

    if tmp_file.stat().st_size <= 1024:
        tmp_file.unlink()
        raise RuntimeError(f"Downloaded file is unexpectedly small: {out_file.name}")
    tmp_file.replace(out_file)


def safe_extract_zip(zip_file: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_root = out_dir.resolve()
    with zipfile.ZipFile(zip_file) as archive:
        for member in archive.infolist():
            target = (out_root / member.filename).resolve()
            if out_root != target and out_root not in target.parents:
                raise RuntimeError(f"Unsafe path inside zip: {member.filename}")
        archive.extractall(out_root)


def ensure_safe_dir(zip_file: Path, out_dir: Path, overwrite: bool = False) -> Path:
    safe_dir = out_dir / zip_file.name.removesuffix(".zip")
    if safe_dir.exists() and not overwrite:
        print(f"Already extracted: {safe_dir}")
        return safe_dir

    with tempfile.TemporaryDirectory(dir=out_dir) as tmp:
        tmp_dir = Path(tmp)
        safe_extract_zip(zip_file, tmp_dir)
        extracted = tmp_dir / safe_dir.name
        if not extracted.exists():
            candidates = list(tmp_dir.glob("*.SAFE"))
            if len(candidates) != 1:
                raise RuntimeError(f"Could not identify extracted SAFE directory in {zip_file}")
            extracted = candidates[0]
        if safe_dir.exists():
            shutil.rmtree(safe_dir)
        shutil.move(str(extracted), safe_dir)
    return safe_dir


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download Sentinel-1 AUX_INS/AUX_CAL files valid for a GAFA RAW product.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--product", help="Sentinel-1 RAW/SLC SAFE path or product name.")
    parser.add_argument("--platform", choices=["S1A", "S1B", "S1C", "S1D"])
    parser.add_argument("--date", help="Acquisition start as YYYYMMDD or YYYYMMDDTHHMMSS.")
    parser.add_argument("--out-dir", default="S1_AUX", help="Directory where AUX files are stored.")
    parser.add_argument(
        "--types",
        nargs="+",
        choices=sorted(AUX_PAGES),
        default=["INS", "CAL"],
        help="AUX types to download. PP1 is only needed for apply_b2b_gain=True.",
    )
    parser.add_argument("--keep-zip", action="store_true", help="Keep downloaded .SAFE.zip files.")
    parser.add_argument("--overwrite", action="store_true", help="Redownload/re-extract existing files.")
    parser.add_argument("--dry-run", action="store_true", help="Print selected files without downloading.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    out_dir = Path(args.out_dir).expanduser().resolve()

    platform, acquisition_time = parse_product_or_date(args.product, args.platform, args.date)
    print(f"Selecting AUX for {platform} at {acquisition_time:%Y%m%dT%H%M%S}")

    for aux_type in args.types:
        products = list_aux_products(aux_type)
        product = choose_aux_product(products, platform, acquisition_time, aux_type)
        zip_file = out_dir / f"{product.name}.SAFE.zip"
        print(
            f"AUX_{aux_type}: {product.name} "
            f"(valid from {product.validity_start:%Y%m%dT%H%M%S})"
        )
        if args.dry_run:
            print(f"  {product.download_url}")
            continue
        download_file(product.download_url, zip_file, overwrite=args.overwrite)
        ensure_safe_dir(zip_file, out_dir, overwrite=args.overwrite)
        if not args.keep_zip:
            zip_file.unlink(missing_ok=True)

    print(f"Done. Set GAFA aux path to: {out_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
