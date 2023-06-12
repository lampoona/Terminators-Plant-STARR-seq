#!/usr/bin/env python

import argparse
import logging
from pathlib import Path
import re
import sys

from pyfaidx import Fasta


G4_REGEX = {
  "AAAA": re.compile(r"(G{3,}).{1,7}\1.{1,7}\1.{1,7}\1|(C{3,}).{1,7}\2.{1,7}\2.{1,7}\2"),
  "AABB": re.compile(r"(G{3,}).{1,7}\1.{0,7}(C{3,}).{1,7}\2"),
  "BBAA": re.compile(r"(C{3,}).{1,7}\1.{0,7}(G{3,}).{1,7}\2"),
  "ABBA": re.compile(r"(G{3,}).{0,7}(C{3,}).{1,7}\2.{0,7}\1|(C{3,}).{0,7}(G{3,}).{1,7}\4.{0,7}\3"),
  "ABAB": re.compile(r"(G{3,}).{0,7}(C{3,}).{0,7}\1.{0,7}\2"),
  "BABA": re.compile(r"(C{3,}).{0,7}(G{3,}).{0,7}\1.{0,7}\2"),
  "ABBB": re.compile(r"(G{3,}).{0,7}(C{3,}).{1,7}\2.{1,7}\2|(G{3,}).{1,7}\3.{1,7}\3.{0,7}(C{3,})"),
  "BAAA": re.compile(r"(C{3,}).{0,7}(G{3,}).{1,7}\2.{1,7}\2|(C{3,}).{1,7}\3.{1,7}\3.{0,7}(G{3,})"),
  "ABAA": re.compile(r"(G{3,}).{0,7}(C{3,}).{0,7}\1.{1,7}\1|(C{3,}).{1,7}\3.{0,7}(G{3,}).{0,7}\3"),
  "BABB": re.compile(r"(C{3,}).{0,7}(G{3,}).{0,7}\1.{1,7}\1|(G{3,}).{1,7}\3.{0,7}(C{3,}).{0,7}\3")
}


def main(genome_path: Path):
  if not genome_path.exists():
      logging.error("%s doesn't exist" % genome_path)
      sys.exit(1)
  genome = Fasta(str(genome_path), as_raw = True)

  for contig in genome:
    for topology_class in G4_REGEX:
      logging.info("Scanning %s for topology class %s" % (contig.name, topology_class))
      for match in G4_REGEX[topology_class].finditer(str(contig)):
        match_start, match_end = match.span()
        sys.stdout.write("\t".join([contig.name, str(match_start), str(match_end), topology_class]) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
      description = "Scan for g-quadruplexes in a genome",
      usage = f"{sys.argv[0]} [options] path/to/genome.fa > output.bed"
    )

    parser.add_argument(
      "genome_path",
      help = "path to genome FASTA to scan",
      type = Path,
      metavar = "path/to/genome.fa"
    )

    parser.add_argument(
      "-v",
      "--verbose",
      help = "output progress and other informative messages",
      action = "count",
      dest = "verbosity",
      default = 0,
    )

    args = parser.parse_args()

    logging.basicConfig(
      level = max(logging.DEBUG, logging.WARNING - (args.verbosity * 10)),
      format = '[{asctime} {levelname}] {message}',
      style = '{'
    )

    main(args.genome_path)
