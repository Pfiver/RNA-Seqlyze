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

    var uploads = {
        'inputfile': {},
        'genbankfile': {}
    };

    var context = function(name) {

        var self = this;

        var options = {
            url:      'upload',
            runtimes:      'html5,gears,flash,silverlight,browserplus,html4',
            browse_button:    name + '_browse',
            drop_element:       name + '_progress',
            multipart_params:     { 'type': name,
                                      'session': upload_session },
            flash_swf_url:          path_js + '/plupload.flash.swf',
            silverlight_xap_url:      path_js + '/plupload.silverlight.xap',
        };

        var events = {
            // 'Init': function(up, params) {
            //     $('#' + name + '_progress .filestatus').text(
            //                     "Current runtime: " + params.runtime);
            // },
            'FilesAdded': function(up, up_files) {
                // remove all other files already present
                // plupload features multiple files in one widget
                // we have two widgets and one name per widget
                up.splice();
                self.active = true;
                $('#' + name + '_progress .filestatus').text(
                    up_files[0].name +
                        ' (' + plupload.formatSize(up_files[0].size) + ')');
            },
            'UploadComplete': function(up, up_files) {
                self.complete = true;
                for (nam in uploads)
                    if (uploads[nam].active)
                        if(!uploads[nam].complete)
                            return;
                $('#create_form').submit();
            },
            'UploadProgress': function(up, up_file) {
                $('#' + name + '_progress .bar').css(
                            "width", up_file.percent + '%');
                // $('#' + name + '_progress .filestatus').text(
                //                             up_file.percent + '%');
            },
        };

        this.active = false;
        this.complete = false;

        var up = this.up = new plupload.Uploader(options);

        $('#' + name + '_progress').click(function() {
            $('#' + name + '_browse').click();
        });

        $('#create_form_submit').click(function() {
            for (nam in uploads)
                if (uploads[nam].active)
                    { up.start(); return false; }
        });

        for (x in events)
            up.bind(x, events[x]);

        up.init();
    };

    for (name in uploads)
        uploads[name] = new context(name);
});

// vim: et:sw=4
