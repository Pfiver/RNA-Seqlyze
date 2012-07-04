" ~/.vimrc - read after /etc/vim/vimrc.local

" let the "sh" syntax plugin know we're using bash
""""""""""""""""""""""""""""""""""""""""""""""""""
let g:is_bash = 1


" highlight over-long lines
"""""""""""""""""""""""""""
highlight OverLength ctermbg=red ctermfg=white
match OverLength /\%81v.\+/

" quit quit quit & stop nagging
"""""""""""""""""""""""""""""""
cabbrev q qa
cabbrev wq wqa
cabbrev wn w<bar>bn
cabbrev wN w<bar>bN
cabbrev n bn
cabbrev N bN

" enable hidden modified buffers
""""""""""""""""""""""""""""""""
set hidden

" quickbuf hotkey
"""""""""""""""""
let g:qb_hotkey = "<F3>"

" language specific settings
""""""""""""""""""""""""""""
function s:python_settings()
	set expandtab
	set tabstop=4
	set shiftwidth=4
endfunction

function s:html_settings()
	set shiftwidth=2
	set indentkeys-=*<Return>
endfunction

function s:javascript_settings()
	set sw=4
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
inoremap <F5> <C-R>=strftime("%a, %d %b %Y %H:%M:%S %z")<CR>

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

autocmd BufEnter,BufNewFile * set indentkeys-=<>>

" simpler code pasting
""""""""""""""""""""""
" insert mode
set pastetoggle=<F2>
" shows new mode when toggling in normal mode
nnoremap <F2> :set invpaste paste?<CR>

" fix auto-indentation
""""""""""""""""""""""
" autocmd BufEnter,BufNewFile * set indentkeys=

" navigate to next/prev file in ":args" using Enter/Backspace
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"nnoremap <silent> <Return> :bn<CR>
"nnoremap <silent> <Backspace> :bN<CR>
nnoremap <silent> <Return> :call g:maybe_next_buffer()<CR>
function g:maybe_next_buffer()
    if bufnr("%") != bufnr('$')
        bn
    endif
endfunction
nnoremap <silent> <Backspace> :call g:maybe_previous_buffer()<CR>
function g:maybe_previous_buffer()
	redir @y | silent ls! | redir END
	if bufnr("%") != matchstr(@y, '[0-9]\+ ') + 0
		bN
	endif
endfunction

" highlight all search pattern matches
""""""""""""""""""""""""""""""""""""""
set hlsearch

" [space] => turn off all highlights
""""""""""""""""""""""""""""""""""""
nnoremap <silent> <Space> :nohlsearch<Bar>:echo<CR>

" resore the cursor position when opening files
"""""""""""""""""""""""""""""""""""""""""""""""
"autocmd BufWinEnter * if ! exists("s:NewFile") | exe "normal g'\"" | endif

" auto-load templates
"""""""""""""""""""""
autocmd BufNewFile * silent! exe "0r ~/.vim/templates/".&ft.".%:e | let s:NewFile=1 | $"

" vim:ts=4:sw=4
