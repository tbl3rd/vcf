(let* ((directory (concat default-directory "src")))
  (mapcar #'cider-load-file (directory-files directory t ".*\\.clj")))
