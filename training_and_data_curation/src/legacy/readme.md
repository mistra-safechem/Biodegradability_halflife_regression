# Legacy modules
Due to refactoring for the intervall based API, modules as used by the p
oint-based invtervall API are stored here.

This way no major refactoring was required, which could have been possible 
but deemed unnecessrary for the time being (time/resources) for something that is not equally interesting.

For some of the earlier "legacy_optimization_scripts", of note:

when using db_utils, ensure you replace db_utils.get_all_data() with db_utils.get_selected_data()
in case it hasn't already been refactored to that.
