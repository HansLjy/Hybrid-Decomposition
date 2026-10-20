\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{weeklyprogress}[2026/03/12 Weekly mathematical progress report]

\LoadClass[11pt]{article}

\RequirePackage[margin=0.75in]{geometry}
\RequirePackage{amsmath}
\RequirePackage{amssymb}
\RequirePackage{amsthm}
\RequirePackage{bm}
\RequirePackage{mathtools}
\RequirePackage{tikz}
\RequirePackage{enumitem}
\RequirePackage{microtype}

\setlength{\parindent}{0pt}
\setlength{\parskip}{0.5em}
\setlist{nosep}

\newtheorem{theorem}{Theorem}[section]
\newtheorem{proposition}[theorem]{Proposition}
\newtheorem{lemma}[theorem]{Lemma}
\newtheorem{corollary}[theorem]{Corollary}
\theoremstyle{definition}
\newtheorem{definition}[theorem]{Definition}
\newtheorem{example}[theorem]{Example}
\newtheorem{algorithm}[theorem]{Algorithm}
\theoremstyle{remark}
\newtheorem*{remark}{Remark}

\newcommand{\meetingdate}[1]{\date{#1}}
\newcommand{\completedwork}{\section{Work Completed}}
\newcommand{\furtherquestions}{\section{Further Questions}}

\AtBeginDocument{\maketitle}
