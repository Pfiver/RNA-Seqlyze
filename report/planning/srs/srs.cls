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

\usepackage[utf8]{inputenc}	% Umlauts
\usepackage[german]{babel}	% ??? Line-breaks

\usepackage{graphicx}		% includegraphics[widdth=...]
\graphicspath{{.}}
\graphicspath{{./img/}}
\usepackage{scrtime}		% ??? nice \today output
\usepackage{microtype}		% optimierter text-satz

\usepackage{xcolor}		% tables
\usepackage{booktabs}
\usepackage{colortbl}
\usepackage{tabularx}

\usepackage{hyperref}		% hyperlinks & -references

\setlength{\parindent}{0pt}	% Keine Abschnitts-Einrückung

\usepackage{scrpage2}		% headers & footers
\usepackage{lastpage}           % needed for "page x of y" headers/footers
\pagestyle{scrheadings}
\clearscrheadfoot
\setheadsepline{1pt}
\renewcommand{\sectionmark}[1]{\markright{#1}{}}
\chead{\rightmark}
\cfoot{\thepage\ / \pageref*{LastPage}}

% font
%\usepackage[T1]{fontenc}
%\usepackage{fourier}
%\usepackage{textcomp}
%\usepackage{lmodern}
%\usepackage[T1]{fontenc}

%\renewcommand{\sfdefault}{pefri8a}
%\renewcommand{\ttdefault}{pef}
\renewcommand{\familydefault}{\sfdefault}


%\usepackage{sectsty}				% Custom sectioning (see below)
%\allsectionsfont{%				% Change font of al section commands
%	\usefont{OT1}{pef}{r}{i}%		% bch-b-n: CharterBT-Bold font
%	\usefont{OT1}{bch}{b}{n}%		% bch-b-n: CharterBT-Bold font
%	\hspace{15pt}%				% Uncomment for indentation
%	}
%\sectionfont{%					% Change font of \section command
%	\usefont{OT1}{bch}{b}{n}%		% bch-b-n: CharterBT-Bold font
%	\sectionrule{0pt}{0pt}{-5pt}{0.8pt}%	% Horizontal rule below section
%	}

\ifdraft
	%\usepackage{lineno}		% Development helpers
	%	\linenumbers		% display line numbers
	%	\modulolinenumbers[2]
	\usepackage{todonotes}		% display todo notes
	%\usepackage{savetrees}		% minimal page margins
\else
	\usepackage[disable]{todonotes}	% hide todo notes
\fi

