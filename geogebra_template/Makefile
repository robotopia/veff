scintarcs.ggb: geogebra_defaults2d.xml geogebra_defaults3d.xml geogebra_javascript.js geogebra.xml images/ss.png
	zip $@ $^

geogebra.xml: geogebra_template.xml values.txt
	cp $< $@
	while read LINE; do \
		TOKEN=`echo $$LINE | awk '{print $$1}'`; \
		VALUE=`echo $$LINE | awk '{print $$2}'`; \
		sed -i "s/$$TOKEN/$$VALUE/g" $@; \
	done < $(word 2,$^)

