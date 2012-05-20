" ~/.vimrc - read after /etc/vim/vimrc.local

" language specific settings
""""""""""""""""""""""""""""
function s:python_settings()
	set expandtab
	set tabstop=4
	set shiftwidth=4
endfunction

" set language specific settings
""""""""""""""""""""""""""""""""
autocmd BufEnter,BufNewFile *
\	if exists("*s:".&ft."_settings")	|
\		exe "call s:".&ft."_settings()"	|
\	endif

" dark background
"""""""""""""""""
set background=dark

" auto-load file specific plugins & indentation heuristics
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
filetype plugin indent on

" F5 -> paste date-time dpkg changelog format
"""""""""""""""""""""""""""""""""""""""""""""
:inoremap <F5> <C-R>=strftime("%a, %d %b %Y %H:%M:%S %z")<CR>

" syntax highlithing
""""""""""""""""""""
syntax on

" various
"""""""""
set modeline " modelines - e.g. // vim:ts=2:sw=2
set smarttab " make <Tab> insert 'shiftwidth' spaces at beginning of line
set incsearch " incremental search
set comments= " ...
set nofoldenable " disable code folding

" fix auto-indentation
""""""""""""""""""""""
" autocmd BufEnter,BufNewFile * set indentkeys=

" navigate to next/prev file in ":args" using Enter/Backspace
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:nnoremap <silent> <Return> :n<CR>
:nnoremap <silent> <Backspace> :N<CR>

" highlight all search pattern matches
""""""""""""""""""""""""""""""""""""""
set hlsearch

" [space] => turn off all highlights
""""""""""""""""""""""""""""""""""""
:nnoremap <silent> <Space> :nohlsearch<Bar>:echo<CR>

" resore the cursor position when opening files
"""""""""""""""""""""""""""""""""""""""""""""""
autocmd BufWinEnter * if ! exists("s:NewFile") | exe "normal g'\"" | endif

" auto-load templates
"""""""""""""""""""""
autocmd BufNewFile * silent! exe "0r ~/.vim/templates/".&ft.".%:e | let s:NewFile=1 | $"

" vim:ts=4:sw=4
