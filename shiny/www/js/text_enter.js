$(document).keyup(function(event) {
    if ($("#type_cpop").is(":focus") && (event.key == "Enter")) {
        $("#go_cpop").click();
    }
});