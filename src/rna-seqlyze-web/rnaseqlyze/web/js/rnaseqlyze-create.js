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

    var uploaders = {
        'inputfile': {},
        'genbankfile': {},
    }

    function upload_complete(up, up_files) {

        uploads[up.id].complete = true;

        complete = true;
        for (name in uploaders)
            if (!uploaders[name].complete) {
                complete = false;
                break;
            }

        if (complete)
            $('#create_form').submit();

    }

    for (name in uploaders) (function(name) {

        options = {
            id:     name,
            url:      'upload',
            runtimes:      'html5,gears,flash,silverlight,browserplus,html4',
            browse_button:    name + '_browse',
            drop_element:       name + '_progress',
            multipart_params:     { 'type': name,
                                      'session': upload_session },
            flash_swf_url:          path_js + '/plupload.flash.swf',
            silverlight_xap_url:      path_js + '/plupload.silverlight.xap',
        };

        var ul = uploaders[name];

        ul.complete = false;

        ul = ul.plupload = new plupload.Uploader(options);

        /*
            ul.bind('Init', function(up, params) {
                $('#' + name + '_container .filestatus').text(
                                "Current runtime: " + params.runtime);
            });
        */

        ul.bind('FilesAdded', function(up, up_files) {
            // remove all other files already present
            // plupload features multiple files in one widget
            // we have two widgets and one name per widget
            up.splice();
            $('#' + name + '_progress .filestatus').text(
                up_files[0].name +
                    ' (' + plupload.formatSize(up_files[0].size) + ')');
        });

        ul.bind('UploadComplete', upload_complete);

        /*    
            ul.bind('UploadProgress', function(up, up_file) {
                $('#' + name + '_progress .bar').css(
                            "width", up_file.percent + '%');
                $('#' + name + '_progress .filestatus').text(
                                            up_file.percent + '%');
            });
        */

        $('#' + name + '_container .progress').click(function() {
            $('#' + name + '_browse').click();
        });
    
        $('#create_form_submit').click(function() {
            console.log(name + " start...")
            uploader.start();
            return false;
        });

        // console.log(name + "_uploader init " + ul.id);
        
        ul.init();

    })(name);
});

// vim: et:sw=4
