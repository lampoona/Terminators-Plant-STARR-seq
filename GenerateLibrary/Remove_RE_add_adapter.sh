awk -v OFS='\t' '{
    NAME=$1
    getline
    SEQ="TAAG"$1"AGGT"
    MUT=""
    if(SEQ ~ /((G|C)GTCTC)|(GAGAC(C|G))/) {
      while(index(SEQ, "GGTCTC") > 0) {MUT=MUT "T" index(SEQ, "GGTCTC") -2 ">A;"; sub(/GGTCTC/, "GGACTC", SEQ)};
      while(index(SEQ, "CGTCTC") > 0) {MUT=MUT "T" index(SEQ, "CGTCTC") - 2 ">A;"; sub(/CGTCTC/, "CGACTC", SEQ)};
      while(index(SEQ, "GAGACC") > 0) {MUT=MUT "A" index(SEQ, "GAGACC") - 1 ">T;"; sub(/GAGACC/, "GAGTCC", SEQ)};
      while(index(SEQ, "GAGACG") > 0) {MUT=MUT "A" index(SEQ, "GAGACG") - 1 ">T;"; sub(/GAGACG/, "GAGTCC", SEQ)};
      print substr(NAME, 2), substr(MUT, 1, length(MUT) - 1) > "RE_mutations_arabidopsis.tsv"
    }
    SEQ="GCGCCGTCTCC"SEQ"CGAGACGGTGC"
    print NAME"\n"SEQ
  }' Arabidopsis_PACs_Thomas_plus_TAIR.fa \
  > Arabidopsis_PACs_Thomas_plus_TAIR_final.fa


awk -v OFS='\t' '{
    NAME=$1
    getline
    SEQ="TAAG"$1"AGGT"
    MUT=""
    if(SEQ ~ /((G|C)GTCTC)|(GAGAC(C|G))/) {
      while(index(SEQ, "GGTCTC") > 0) {MUT=MUT "T" index(SEQ, "GGTCTC") -2 ">A;"; sub(/GGTCTC/, "GGACTC", SEQ)};
      while(index(SEQ, "CGTCTC") > 0) {MUT=MUT "T" index(SEQ, "CGTCTC") - 2 ">A;"; sub(/CGTCTC/, "CGACTC", SEQ)};
      while(index(SEQ, "GAGACC") > 0) {MUT=MUT "A" index(SEQ, "GAGACC") - 1 ">T;"; sub(/GAGACC/, "GAGTCC", SEQ)};
      while(index(SEQ, "GAGACG") > 0) {MUT=MUT "A" index(SEQ, "GAGACG") - 1 ">T;"; sub(/GAGACG/, "GAGTCC", SEQ)};
      print substr(NAME, 2), substr(MUT, 1, length(MUT) - 1) > "RE_mutations_maize.tsv"
    }
    SEQ="GCGCCGTCTCC"SEQ"CGAGACGGTGC"
    print NAME"\n"SEQ
  }' Maize_Top_PACs_Final.fa \
  > Maize_Top_PACs_Final_Final.fa
