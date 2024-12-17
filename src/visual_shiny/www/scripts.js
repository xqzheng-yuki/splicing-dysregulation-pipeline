$(function () {
    $('[data-toggle=tooltip]').tooltip()
  })
  
$(document).keyup(function(event) {
    if ($("#gene").is(":focus") && (event.key == "Enter")) {
        $("#get_gene").click();
    }
});