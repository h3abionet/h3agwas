head -n 8 $manifest > m.csv
grep $1 $manifest >> m.csv
head -n 6 $strand   > s.csv
grep $1 $strand >> s.csv