## Makefile di esempio per il progetto di High Performance Computing
## 2019/2020, corso di laurea in Ingegneria e Scienze Informatiche,
## Universita' di Bologna.

## Ultima modifica: 2019-11-16, Moreno Marzolla <moreno.marzolla(at)unibo.it>

## Questo e' un frammento del Makefile utilizzato negli esempi
## illustrati durante il corso. Questo file puo' essere modificato e
## personalizzato in base alle proprie esigenze. Si puo' anche
## decidere di non usarlo; in tal caso indicare le istruzioni di
## compilazione nel file README presente nella directory a livello
## superiore.
##
## Se si decide di usare questo makefile, il comando "make" dovrebbe
## compilare tutti i programmi consegnati senza errori né warning.  Si
## consiglia pertanto di rimuovere eventuali target non presenti
## nell'archivio consegnato.
##
## Questo Makefile compila i file "omp-*.c" usando il flag -fopenmp, i
## file "cuda-*.cu" con il compilatore nvcc, e i file "mpi-*.c" con
## mpicc.
##
## I principali target definiti da questo makefile sono:
##
## make         compila tutti i sorgenti disponibili
## make clean   cancella i file temporanei e gli eseguibili
## make openmp  compila la versione OpenMP
## make mpi     compila la versione MPI
## make cuda    compila la versione CUDA
## make images  usa ./convex-hull per creare le immagini di tutti gli inviluppi

SRC_DIR:=src/
OUT_DIR:=build/
INPUTS_DIR:=inputs/
IMAGES_DIR:=images/
DOCS_DIR:=docs/

EXE_OMP:=$(addprefix ${OUT_DIR}, $(basename $(notdir $(wildcard ${SRC_DIR}omp-*.c))))
EXE_MPI:=$(addprefix ${OUT_DIR}, $(basename $(notdir $(wildcard ${SRC_DIR}mpi-*.c))))
EXE_CUDA:=$(addprefix ${OUT_DIR}, $(basename $(notdir $(wildcard ${SRC_DIR}cuda-*.cu))))
EXE_SERIAL:=$(addprefix ${OUT_DIR}, convex-hull)
DATAFILES:=$(addprefix ${INPUTS_DIR}, ace.in box1k.in box10k.in box100k.in box1M.in circ1k.in circ10k.in circ100k.in gaus100k.in gaus1M.in)
IMAGES:=$(patsubst ${INPUTS_DIR}%.in, ${IMAGES_DIR}%.png, $(DATAFILES))
EXE:=$(EXE_OMP) $(EXE_MPI) $(EXE_SERIAL) $(EXE_CUDA)
CFLAGS+=-std=c99 -Wall -Wpedantic -O1 -D_XOPEN_SOURCE=600
LDLIBS+=-lm
NVCC:=nvcc
NVCFLAGS+=-Wno-deprecated-gpu-targets
NVLDLIBS+=-lm
REPORT:=$(DOCS_DIR)report.pdf

.PHONY: clean

ALL: $(EXE) datafiles

datafiles: $(INPUTS_DIR) $(DATAFILES)

images: $(IMAGES_DIR) $(IMAGES)

docs: $(REPORT)

% : %.cu
	$(NVCC) $(NVCFLAGS) $< -o $@ $(NVLDLIBS)

$(OUT_DIR)%: $(SRC_DIR)%.c
	\mkdir -p $(OUT_DIR)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $< -o $@ $(LOADLIBES) $(LDLIBS)

$(EXE_OMP): CFLAGS+=-fopenmp
$(EXE_OMP): LDLIBS+=-lgomp
openmp: $(EXE_OMP)

$(EXE_MPI): CC=mpicc
mpi: $(EXE_MPI)

cuda: $(EXE_CUDA)

$(INPUTS_DIR)box1k.in:
	rbox 1000 D2 > $@

$(INPUTS_DIR)box10k.in:
	rbox 10000 D2 > $@

$(INPUTS_DIR)box100k.in:
	rbox 100000 D2 > $@

$(INPUTS_DIR)box1M.in:
	rbox 1000000 D2 > $@

$(INPUTS_DIR)box10M.in:
	rbox 10000000 D2 > $@

$(INPUTS_DIR)circ1k.in:
	rbox s 1000 D2 > $@

$(INPUTS_DIR)circ10k.in:
	rbox s 10000 D2 > $@

$(INPUTS_DIR)circ100k.in:
	rbox s 100000 D2 > $@

%.hull: %.in $(EXE_SERIAL)
	echo create $@
	./$(EXE_SERIAL) < $< > $@

$(IMAGES_DIR)%.png: $(INPUTS_DIR)%.in $(INPUTS_DIR)%.hull
	gnuplot -c plot-hull.gp $+ $@

$(DOCS_DIR)%.pdf: $(filter-out %.log %.toc %.aux %.pdf, $(wildcard $(DOCS_DIR)*))
	pdflatex -halt-on-error -output-directory $(DOCS_DIR) $(@:.pdf=.tex)

%/:
	mkdir $@

clean:
	\rm -rf $(OUT_DIR) *.o *~ $(IMAGES_DIR) $(INPUTS_DIR)*.hull
