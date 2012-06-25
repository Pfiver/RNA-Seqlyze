/*
 * RNA-seqlyze javascript routines
 *
 *     for the "create" page
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


    /*
     * plupload -- from src/plupload/examples/custom.html
     */
    
    jQuery.each(['inputfile', 'genbankfile'], function(i, file) { (function(file) {
        
        var uploader = new plupload.Uploader({
            runtimes:           'html5,gears,flash,silverlight,browserplus,html4',
            browse_button:      file + '_browse',
            drop_element:       file + '_progress',
            url:                'upload',
            multipart_params:   { 'file': file },
            flash_swf_url:      path_js + '/plupload.flash.swf',
            silverlight_xap_url:    path_js + '/plupload.silverlight.xap',
        });

        /*
            uploader.bind('Init', function(up, params) {
                $('#' + file + '_container .filestatus').text("Current runtime: " + params.runtime);
            });
        */

        uploader.bind('FilesAdded', function(up, files) {
            up.splice();    // removes all other files already present
            $('#' + file + '_progress .filestatus').text(files[0].name + ' (' + plupload.formatSize(files[0].size) + ')');
        });
    
        uploader.bind('UploadProgress', function(up, file) {
            $('#' + file + '_progress .filestatus').text(file.percent + '%');
        });
        
        $('#' + file + '_container .progress').click(function() {
            $('#' + file + '_browse').click();
        });
    
        $('#create_form_submit').click(function() {
            console.log(file + " start...")
            uploader.start();
            return false;
        });

        // console.log(file + "_uploader init " + uploader.id);
        
        uploader.init();
        
    })(file);});

});
