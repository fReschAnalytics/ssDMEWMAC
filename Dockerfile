FROM continuumio/anaconda3

ADD ssDMEWMAC/__init__.py /ssDMEWMAC/
ADD ssDMEWMAC/functions.py /ssDMEWMAC/
ADD data/h_init_lookup.csv /data/
ADD data/h_lookup.csv /data/
ADD in_control_simulation.py /
ADD out_of_control_simulation.py /
ADD control.sh /

RUN pip install numpy pandas scipy scikit-learn gcloud

RUN mkdir /out

CMD [ "bash", "control.sh"]