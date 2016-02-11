all:
	lualatex writeup
	bibtex writeup
	lualatex writeup
	lualatex writeup
clean:
	-@rm -rf *.log
	-@rm -rf *.aux
	-@rm -rf writeupNotes.bib
	-@rm -rf writeup.bbl
	-@rm -rf writeup.toc
	-@rm -rf writeup.out
	-@rm -rf writeup.blg
	-@rm -rf writeup.nav
	-@rm -rf writeup.snm

open:
	gnome-open writeup.pdf </dev/null &>/dev/null &
# TODO get the correct viewer for the laptop installed, check
