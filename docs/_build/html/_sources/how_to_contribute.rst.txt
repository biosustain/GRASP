How to contribute
=======================

Contributions are always welcome! :)

Below are some guidelines you should follow to make sure no bugs are introduced in the process (or at least are minimized). 


Get started!
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here's the workflow you should follow when implementing any changes in the code and/or documentation:

1. If you are not a collaborator in the repo, start by forking it.

2. Clone either your fork (if you are **not** a collaborator on the repo) or the main repo (if you are a collaborator on the repo) into your machine.
  
   .. code-block:: console
        
        git clone <repo_url>


3. Create a local branch starting from `master`, where you will implement your changes. Remember to give it a meaningful name.
  
   .. code-block:: console

        git checkout master
        git branch <branch_name>
        git checkout <branch_name>


4. When making your changes commit often. This will make it easier to review your work. 
Please use `semantic commit messages <https://karma-runner.github.io/2.0/dev/git-commit-msg.html>`_.

   .. code-block:: console

        git add <file_name1>
        git commit


5. Remember to either update the unit tests for any functions that you changed or add unit tests for any functions that you added.

6. Run the unit tests by going to ``matlab_code/tests`` and running ``run_unit_tests``.

7. Update or add relevant documentation.

8. If all the unit tests pass, go ahead and push your branch to Github:
   
   .. code-block:: console

        git push origin <branch_name>


9. On Github, create a pull request from your branch to ``master``, provide a detailed description of the changes introduced, and wait for feedback.


Please try to follow the existing coding style.


Reporting bugs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To report a bug, please create an issue at https://github.com/biosustain/GRASP/issues.






