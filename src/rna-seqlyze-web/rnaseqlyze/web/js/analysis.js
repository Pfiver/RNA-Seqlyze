/*
 * RNA-seqlyze javascript routines
 */

$(document).ready(function() {

    /*
     * function library
     */

    function maybe_show_pairendlen_controls() {
        if ($('#pairendedInput').attr('checked'))
            $('#pairendlenControls').show();
        else
            $('#pairendlenControls').hide();
    }

    function search_organism() {
        oc = $('#organismControls');
        os = $('#organismStatus');
        oi = $('#organismInput');

        oc.removeClass("success");
        oc.removeClass("error");

        if (oi.val().match(/^NC/)) {
            oc.addClass("success");
            os.text("found");
        } else {
            oc.addClass("error");
            os.text("not found");
        }
    }

    /*
     * page initialization
     */

    maybe_show_pairendlen_controls();
    $('#pairendedInput').change(maybe_show_pairendlen_controls);

    $('#organismInput').blur(search_organism);

});
