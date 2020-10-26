f2py=f2py
python=python

all: TMscoreMod nwalign ScoreMod tmalign

TMscoreMod: PDB/TMscore/TMscoreMod.f
	cd PDB/TMscore;${f2py} -c -m TMscoreMod TMscoreMod.f

nwalign: PDB/NWalign/NWalign_wrapper.f
	cd PDB/NWalign;${f2py} -c -m nwalign NWalign_wrapper.f

ScoreMod: PDB/KABSCH_FORTRAN/ScoreMod.f
	cd PDB/KABSCH_FORTRAN;${f2py} -c -m ScoreMod ScoreMod.f

tmalign: PDB/TMalign/TMalign_wrapper.h PDB/TMalign/TMalign_wrapper.cpp
	cd PDB/TMalign;python setup.py build_ext --inplace

clean:
	rm PDB/*/*.so
