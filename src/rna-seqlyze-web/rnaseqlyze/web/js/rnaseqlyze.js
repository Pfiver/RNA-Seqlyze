/*
 * RNA-seqlyze javascript routines
 */

$(function() {

    // global variables and helpers
    // ----------------------------

    // ___ el ___
    //  from http://joestelmach.github.com/laconic/
    window.el = $.el;

    // log.info() and log.debug()
    window.log = {
        'info': console.log,
        'debug': function() {}
    }
    if (rnaseqlyze_debug)
        window.log.debug = console.log;

    // page initializaion
    // ------------------

    // use bootstrap's
    // "scrollspy" plugin
    //  -- patched version - see https://github.com/twitter/bootstrap/pull/3829
    $(window).scrollspy({
//            offset: 200,
            wrap: $('#wrap')[0],
    });

});
