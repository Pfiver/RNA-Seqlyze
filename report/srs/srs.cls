% provide the "srs" documentclass class
\ProvidesClass{srs}[2012/05/23 Software Requirement Specification]

% declare some options
\newif\ifdraft
\draftfalse
\DeclareOption{draftmode}{%
	\drafttrue%
}

% pass unknow options on to base class
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{scrartcl}}

% important - don't move :-)
\ProcessOptions

% This style is based on KOMA-Script
\LoadClass{scrartcl}

% ... but, we're not printing a book, so make bottom=top-margin
\makeatletter
	\advance\textheight by \ta@vblk
	\advance\textheight by -\headsep
\makeatother

% ! never blindy change the order of package imports !

\usepackage{graphicx}		% includegraphics[widdth=...]
\graphicspath{{./}}
\graphicspath{{../../assets/}}
\usepackage{scrtime}		% ??? nice \today output

\usepackage{xcolor}		% tables
\usepackage{booktabs}
\usepackage{colortbl}
\usepackage{tabularx}

\usepackage{hyperref}		% hyperlinks & -references
\hypersetup{
	colorlinks,%
	urlcolor=black,%
	citecolor=black,%
	filecolor=black,%
	linkcolor=black%
}
% in a current version of the package,
% the above lines could be replaced with
%\usepackage[hidelinks]{hyperref}

%\setlength{\parindent}{0pt}	% No paragraph indentaion

\usepackage{scrpage2}		% headers & footers
\usepackage{lastpage}           % needed for "page x of y" headers/footers
\pagestyle{scrheadings}
\clearscrheadfoot
\setheadsepline{1pt}
\renewcommand{\sectionmark}[1]{\markright{\thesection\ - #1}{}}
\chead{\rightmark}
\cfoot{\thepage\ / \pageref*{LastPage}}
\ifoot[\includegraphics{fhnw_hls_e_10mm}]{}

\usepackage{ifxetex}

\ifxetex
	% FHNW Fonts
	% XeTeX needed -- this doesn't work with pdflatex
	% if it doesn't work with xelatex,
	% try comying the font files to ~/.fonts
	% (mkdir ~/.fonts if that directory doesn't exist)
	\usepackage{fontspec}
	\setsansfont[ExternalLocation=../../assets/,
		BoldFont      =*-Bold,
		ItalicFont    =*-Obl,
		BoldItalicFont=*-BoldObl
		]{UniversLTStd}
	\setromanfont[ExternalLocation=../../assets/,
		UprightFont   =*-Roman,
		BoldFont      =*-Bold,
		ItalicFont    =*-Italic,
		BoldItalicFont=*-Italic
		]{EgyptienneFLTStd}
\else
	\usepackage[utf8]{inputenc}
	\usepackage[T1]{fontenc}
	\usepackage{charter}
\fi

% Scrict text justification
\tolerance=10000
\newcommand\nohyphenate{\linebreak[1]\hbox}

\ifdraft
	%\usepackage{lineno}		% Development helpers
	%	\linenumbers		% display line numbers
	%	\modulolinenumbers[2]
	\usepackage{todonotes}		% display todo notes
	%\usepackage{savetrees}		% minimal page margins
\else
	\usepackage[disable]{todonotes}	% hide todo notes
\fi
