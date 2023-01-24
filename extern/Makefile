PUBDIR := public
ROOTDIR := ..

SOURCES := $(shell find ${ROOTDIR}/src -name "*.F*") $(shell find ${ROOTDIR}/unit -name "*.F*") $(shell find ${ROOTDIR}/user -name "*.F*")
SOURCES := $(filter-out ${ROOTDIR}/src/algorithms/gca/%,$(SOURCES))
PUBLIC := $(subst ${ROOTDIR}/, ${PUBDIR}/, $(SOURCES))

public: $(PUBLIC)
	@cp -r ${ROOTDIR}/tristanVis/ ${PUBDIR}/
	@cp -r ${ROOTDIR}/.vscode ${PUBDIR}/
	@cp -r ${ROOTDIR}/docker/ ${PUBDIR}/
	@mkdir -p ${PUBDIR}/extern/
	@cp ${ROOTDIR}/extern/makedepf90 ${PUBDIR}/extern/
	@cp -r ${ROOTDIR}/inputs/ ${PUBDIR}/
	@cp ${ROOTDIR}/Makefile.in ${PUBDIR}/
	@cp ${ROOTDIR}/README.md ${PUBDIR}/
	@cp ${ROOTDIR}/dev-requirements.txt ${PUBDIR}/
	@cp ${ROOTDIR}/fulltest.py ${PUBDIR}/
	@cp ${ROOTDIR}/configure.py ${PUBDIR}/
	@cp ${ROOTDIR}/.gitignore ${PUBDIR}/
	@cp ${ROOTDIR}/LICENSE ${PUBDIR}/
	@cp ${ROOTDIR}/CITATION ${PUBDIR}/

${PUBDIR}/%: ${ROOTDIR}/%
	@mkdir -p $(dir ${@})
	@-./unifdef $< -o $@ -U QED -U BWPAIRPRODUCTION -U COMPTONSCATTERING -U PAIRANNIHILATION -U GCA -B -t

clean:
	@echo "Cleaning ..."
	@rm -rf ${PUBDIR}
