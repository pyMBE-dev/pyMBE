#!/bin/bash

FILES = $(shell grep -Po "^logo_[^:]+(?=:)" Makefile)

all:
	make $(FILES)

list:
	@echo $(FILES)

logo_plain.svg: logo.svg Makefile
	inkscape --actions="export-filename:$@; export-area-page; export-text-to-path; export-plain-svg; export-do;" $<

logo_sepia.png: logo.svg Makefile
	inkscape --actions="export-filename:$@; export-area-page; export-dpi:192; export-do;" $<

logo_sepia.eps: logo.svg Makefile
	inkscape --actions="export-filename:$@; export-area-page; export-do;" $<

logo_sepia.pdf: logo.svg Makefile
	inkscape --actions="export-filename:$@; export-area-page; export-text-to-path; export-do;" $<

logo_white.png: logo.svg Makefile
	inkscape --actions="export-filename:$@; export-area-page; export-id:layer1; export-id-only; export-background:white; export-dpi:192; export-do;" $<

logo_white.eps: logo.svg Makefile
	inkscape --actions="export-filename:$@; export-area-page; export-id:layer1; export-id-only; export-background:white; export-do;" $<

logo_white.pdf: logo.svg Makefile
	inkscape --actions="export-filename:$@; export-area-page; export-id:layer1; export-id-only; export-background:white; export-do;" $<

logo_transparent.png: logo.svg Makefile
	inkscape --actions="export-filename:$@; export-area-page; export-id:layer1; export-id-only; export-dpi:192; export-do;" $<

temp_logo_banner_base.svg: logo.svg Makefile
	inkscape --actions="export-filename:$@; select-by-id:text286611; transform-translate:700,-205; export-id:layer1; export-id-only; export-plain-svg; export-do;" $<
	sed -i '/font-size:22.5778px/,$${s//font-size:67.60740px/;b};$$q1' $@

logo_banner.png: temp_logo_banner_base.svg Makefile
	inkscape --actions="export-filename:$@; export-id:layer1; export-dpi:48; export-do;" $<

logo_banner_dark_mode.png: temp_logo_banner_base.svg Makefile
	cp $< $@.svg
	sed -i '/text-anchor:middle;fill:#3d2c25/,$${s//text-anchor:middle;fill:#eaeaea/;b};$$q1' $@.svg
	inkscape --actions="export-filename:$@; export-id:layer1; export-dpi:48; export-do;" $@.svg
	rm $@.svg
