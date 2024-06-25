PUBDIR := public
ROOTDIR := ..

SOURCES := $(shell find ${ROOTDIR}/src -name "*.F*") $(shell find ${ROOTDIR}/unit -name "*.F*") $(shell find ${ROOTDIR}/user -name "*.F*")
SOURCES := $(filter-out ${ROOTDIR}/src/algorithms/gca/%,$(SOURCES))
PUBLIC := $(subst ${ROOTDIR}/, ${PUBDIR}/, $(SOURCES)) ${PUBDIR}/CMakeLists.txt ${PUBDIR}/src/CMakeLists.txt

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

${PUBDIR}/CMakeLists.txt: ${ROOTDIR}/CMakeLists.txt
	@cp $< $@
	@sed -i.bak '/BW pair production/d' $@
	@sed -i.bak '/QED step/d' $@
	@sed -i.bak '/Compton scattering/d' $@
	@sed -i.bak '/Pair annihilation/d' $@
	@sed -i.bak '/Enable QED/d' $@
	@sed -i.bak '/Enable Breit-Wheeler pair production/d' $@
	@sed -i.bak '/Enable pair annihilation/d' $@
	@sed -i.bak '/Number of GCA/d' $@

${PUBDIR}/src/CMakeLists.txt: ${ROOTDIR}/src/CMakeLists.txt
	@cp $< $@
	@sed -i.bak '/gca/d' $@

clean:
	@echo "Cleaning ..."
	@rm -rf ${PUBDIR}
