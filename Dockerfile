FROM jupyter/scipy-notebook

RUN pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org jupyter_contrib_nbextensions
RUN jupyter contrib nbextension install --user
RUN jupyter nbextension enable spellchecker/main