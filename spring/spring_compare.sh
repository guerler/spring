#!/bin/bash

./scripts/performance/verified_concat.sh /home/guerler/package/results/spring_refinement/features/ results/spring.txt
./scripts/performance/verified_concat.sh /home/guerler/package/results/spring_tasser/features/ results/spring_tasser.txt
./scripts/performance/verified_concat.sh /home/guerler/package/results/spring_classic/features/ results/spring_classic.txt

# compare
./scripts/performance/compare.sh results/spring.txt results/spring_tasser.txt > results/compare_tasser.txt
./scripts/performance/compare.sh results/spring.txt results/spring_classic.txt > results/compare_classic.txt

#./scripts/performance/spring_coverage.sh results/spring.txt > results/spring.coverage.txt
#./scripts/performance/spring_ranking_rmsd.sh results/spring.txt > results/spring.ranking.rmsd.txt
#./scripts/performance/spring_ranking_tmscore.sh results/spring.txt > results/spring.ranking.tmscore.txt
#./scripts/performance/spring_ranking_native.sh results/spring.txt > results/spring.ranking.native.txt
#./scripts/performance/spring_ranking_globalrmsd.sh results/spring.txt > results/spring.ranking.globalrmsd.txt