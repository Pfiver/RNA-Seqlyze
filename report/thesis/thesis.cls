% provide the "thesis" documentclass class
\ProvidesClass{thesis}[2012/07/09 Bachelor Thesis]

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
\graphicspath{
	{./}
	{./img/}
	{../assets/}
}
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
\renewcommand{\sectionmark}[1]{\markleft{\thesection\ - #1}{}}
\def\empty{}\def\trigger{}\def\nextmark{}\def\curmark{}
\renewcommand{\paragraphmark}[1]{\def\nextmark{#1}}
\renewcommand\rightmark{
	\ifx\trigger\empty
	\ifx\nextmark\empty
	\else
		\global\def\trigger{x}
	\fi
	\else
		\scriptsize\curmark
		\global\let\curmark\nextmark
	\fi
}
\ihead{\leftmark}
\ohead{\rightmark}
\cfoot{\thepage\ / \pageref*{LastPage}}

\usepackage{ifxetex}

\ifxetex
	\ifoot[
		\upshape \setlength\parindent{0pt}
		\hspace{-0.34cm} \includegraphics{FHNW_HLS_10mm} \par
		{ \sffamily \scriptsize Institut für Chemie und Bioanalytik \par }
		\vskip -\baselineskip \hfill Muttenz, Juli 2012
	\vspace{999pt}
	]{}

	% FHNW Fonts
	% XeTeX needed -- this doesn't work with pdflatex
	% if it doesn't work with xelatex,
	% try comying the font files to ~/.fonts
	% (mkdir ~/.fonts if that directory doesn't exist)
	\usepackage{fontspec}
	\setsansfont[ExternalLocation=../assets/,
		BoldFont      =*-Bold,
		ItalicFont    =*-Obl,
		BoldItalicFont=*-BoldObl
		]{UniversLTStd}
	\setromanfont[ExternalLocation=../assets/,
		UprightFont   =*-Roman,
		BoldFont      =*-Bold,
		ItalicFont    =*-Italic,
		BoldItalicFont=*-Italic
		]{EgyptienneFLTStd}

	\usepackage[style=alphabetic]{biblatex}

	\addbibresource{Genomics}
\else
	\usepackage[utf8]{inputenc}
	\usepackage[T1]{fontenc}
	\usepackage{textcomp}
	\usepackage{charter}
	\usepackage[ngerman]{babel}

	\ifoot[
		{[FHNW HLS Logo]}
		{ \sffamily \scriptsize Institut für Chemie und Bioanalytik \par }
		\vskip -\baselineskip \hfill Muttenz, Juli 2012
	\vspace{999pt}
	]{}
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
%	\usepackage[disable]{todonotes}	% hide todo notes
	\usepackage{todonotes}		% show todo notes
\fi

\usepackage{longtable}
\usepackage{sphinx}

\renewcommand\figurename{Abb.}

\usepackage[nottoc]{tocbibind}

