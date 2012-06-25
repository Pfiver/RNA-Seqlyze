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

    var uploads = ['inputfile', 'genbankfile'];

    var context = function(name) {
        var self = this;
        this.options = {
            url:      'upload',
            runtimes:      'html5,gears,flash,silverlight,browserplus,html4',
            browse_button:    name + '_browse',
            drop_element:       name + '_progress',
            multipart_params:     { 'type': name,
                                      'session': upload_session },
            flash_swf_url:          path_js + '/plupload.flash.swf',
            silverlight_xap_url:      path_js + '/plupload.silverlight.xap',
        };
        this.events = {
            'Init': function(up, params) {
console.log("1111111 " + name);
                $('#' + name + '_progress .filestatus').text(
                                "Current runtime: " + params.runtime);
            },
            'FilesAdded': function(up, up_files) {
console.log("4444444 " + self.up.id);
console.log("4444445 " + self.up.files);
console.log("4444446 " + up_files);
                // remove all other files already present
                // plupload features multiple files in one widget
                // we have two widgets and one name per widget
                up.splice();
                this.active = true;
                $('#' + name + '_progress .filestatus').text(
                    up_files[0].name +
                        ' (' + plupload.formatSize(up_files[0].size) + ')');
            },
            'UploadComplete': function(up, up_files) {
                this.complete = true;
                for (nam in uploads)
                    if (uploads[nam].active)
                        if(!uploads[nam].complete)
                            return;
                console.log("go");
                $('#create_form').submit();
            },
            'UploadProgress': function(up, up_file) {
                $('#' + name + '_progress .bar').css(
                            "width", up_file.percent + '%');
                // $('#' + name + '_progress .filestatus').text(
                //                             up_file.percent + '%');
            },
        };
        this.init = function() {
            this.up = new plupload.Uploader(this.options);
            for (x in this.events) {
console.log("2222222 " + x);
                this.up.bind(x, this.events[x]);}
            this.up.init();
console.log("5555555 " + self.up.id);
console.log("6666666 " + this.up.id);

            $('#create_form_submit').click(function() {
                console.log(name + " start...")
                self.up.start();
console.log("3333333 " + self.up.id);
                return false;
            });
            $('#' + name + '_progress').click(function() {
                $('#' + name + '_browse').click();
            });
        };
        this.active = false;
        this.complete = false;
    };

    for (i = 0; i < uploads.length; i++) {
        ctx = new context(uploads[i]);
        ctx.init();
    }

});

// vim: et:sw=4
