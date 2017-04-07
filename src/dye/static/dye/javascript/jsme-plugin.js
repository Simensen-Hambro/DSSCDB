//this function will be called after the JavaScriptApplet code has been loaded.
function jsmeOnLoad() {
    $.jsmeApplet = new JSApplet.JSME("jsme_container", "100%", "500px");

    if ($.loaded_smiles) {
        $.jsmeApplet.readGenericMolecularInput($.loaded_smiles);
    }
}

$(function () {
    $.loaded_smiles = $("#id_smiles").val();

    $(".submit-button").click(function (e) {
        e.preventDefault();
        var smiles = $.jsmeApplet.smiles();
        $("#id_smiles").val(smiles);
        console.log(smiles);
        $("#search-form").submit();
    });

    var updateWidth = function () {
        $("#hideFrame").width($("#jsme_container").children().eq(1).width());
        if ($.jsmeApplet) {
            $.jsmeApplet.setSize("100%", "500px");
        }
    };

    setInterval(updateWidth, 200);

    $("#hideFrame").click(function () {
        $("#hideFrame").css("display", "none")
    });
});

