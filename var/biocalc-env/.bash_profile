# tmux
######

tmux="$HOME/.local/bin/tmux -S /tmp/tmux-$USER"

if ! [ $TMUX ]
then
	$tmux display 2> /dev/null && TERM=xterm-256color exec $tmux attach
	TERM=xterm-256color exec $tmux
fi

# other stuff
#############

set -PC

CDPATH=:~

xpath() {
        [[ ${!1} =~ (^|:)$2(:|$) ]] || eval ${1}=\$2:\$$1
}
xpath PATH ~/.local/bin

alias vi=vim

alias cp='cp -i'
alias mv='mv -i'
alias rm='rm -i'
alias grep='grep --color=auto'
alias ls='ls --color=auto --show-control-chars --group-directories-first'

alias gd='git diff'
alias gs='git status'

eval "$(lesspipe)"
eval "$(dircolors -b)"

stty -ixon -ixoff

export VISUAL=vim EDITOR=vim PYTHONSTARTUP=~/.pythonstartup.py

PS1='\[\033[;92m\]\u@\h\[\033[00m\]:\[\033[;94m\]\w\[\033[00m\]\$ '
PROMPT_COMMAND='echo -ne "\033]0;${USER}@$(hostname -s): ${PWD/$HOME/~}\007"'

find-py()
{
	find "$@" \( -name junk -o -name setup.py \) -prune -o \
		\( -name "*.pt" -o \( -name "*.py" \( ! -name __init__.py \
		-o -size +1 -exec sed -n '1{/^__/{q1}}' {} \; \) \) \) -print
}

# viall
#######
# edit all python, javascript, page-templte and lessccss files
# in the RNA-Seqlyze source tree
viall ()
{
    cd $rnas_topdir/src

    not_starts_with__="-exec sed -n 1{/^__/{q1}} {} ;"

    vi $(find rna-seqlyze* \(                           \
            -name junk -o -name setup.py                \
        \) -prune -o \(                                 \
            -name "rnaseqlyze*.js"                      \
            -o -name "rnaseqlyze*.less"                 \
            -o -name "*.pt"                             \
            -o -name "*.py" \(                          \
                -not -name __init__.py                  \
                -o -size +1 $not_starts_with__ \)       \
        \) -print |
                        tr '/-' '[]' | sort | tr '[]' '/-'
    ) 
}

# rna-seqlyze
#############

bt() {
	bt=$rnas_topdir
	wd=$rnas_workdir

	. ~/data/rna-seqlyze/bash-env

	cd "$rnas_topdir"
}
