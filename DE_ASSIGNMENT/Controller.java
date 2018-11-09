package sample;

import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.scene.Node;
import javafx.scene.chart.CategoryAxis;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.scene.chart.XYChart.Data;
import javafx.scene.control.*;
import javafx.scene.input.MouseEvent;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.Pane;
import javafx.scene.shape.MoveTo;
import javafx.scene.shape.Path;
import javafx.stage.Stage;

import java.net.URL;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.ResourceBundle;
import java.util.Set;


public class Controller implements Initializable {
    public Label H;
    @FXML
    private LineChart<Number, Number> chart;
    @FXML
    GridPane gridPane;
    @FXML
    Pane pane;
    @FXML
    NumberAxis xAxis, yAxis;
    @FXML
    Button quitButton, plotButton;
    @FXML
    TextField NInput, NsInput, NfInput, XInput, x0Input, y0Input;
    @FXML
    Label N, Ns, Nf, ch;
    @FXML
    ChoiceBox<String> choiceBox;

    private XYChart.Series<Number, Number> sinSeries, cosSeries;

    public void sayH(ActionEvent act) {
        H.setText("H");
    }

    public void initialize(URL url, ResourceBundle rb) {
        NfInput.setVisible(false);
        NsInput.setVisible(false);
        Ns.setVisible(false);
        Nf.setVisible(false);

        choiceBox.getItems().clear();
        choiceBox.setValue("");
        choiceBox.getItems().addAll("Exact", "Euler", "Euler local error", "Euler total error", "Euler improved","Euler improved local error", "Euler improved total error", "Runge-Kutta", "Runge-Kutta local error", "Runge-Kutta total error");
        choiceBox.getSelectionModel().selectedIndexProperty().addListener(new ChangeListener<Number>() {
            @Override
            public void changed(ObservableValue<? extends Number> observableValue, Number number, Number number2) {
                ch.setVisible(false);
                switch (number2.intValue()) {
                    case 0:
                        NInput.setVisible(true);
                        N.setVisible(true);
                        NInput.setEditable(true);
                        NsInput.setEditable(false);
                        NsInput.setVisible(false);
                        Ns.setVisible(false);
                        NfInput.setEditable(false);
                        NfInput.setVisible(false);
                        Nf.setVisible(false);
                        break;
                    case 1:
                        NInput.setVisible(true);
                        N.setVisible(true);
                        NInput.setEditable(true);
                        NsInput.setEditable(false);
                        NsInput.setVisible(false);
                        Ns.setVisible(false);
                        NfInput.setEditable(false);
                        NfInput.setVisible(false);
                        Nf.setVisible(false);
                        break;
                    case 2:
                        NInput.setVisible(true);
                        N.setVisible(true);
                        NInput.setEditable(true);
                        NsInput.setEditable(false);
                        NsInput.setVisible(false);
                        Ns.setVisible(false);
                        NfInput.setEditable(false);
                        NfInput.setVisible(false);
                        Nf.setVisible(false);
                        break;
                    case 3:
                        NInput.setVisible(false);
                        N.setVisible(false);
                        NInput.setEditable(false);
                        NsInput.setEditable(true);
                        NsInput.setVisible(true);
                        Ns.setVisible(true);
                        NfInput.setEditable(true);
                        NfInput.setVisible(true);
                        Nf.setVisible(true);
                        break;
                    case 4:
                        NInput.setVisible(true);
                        N.setVisible(true);
                        NInput.setEditable(true);
                        NsInput.setEditable(false);
                        NsInput.setVisible(false);
                        Ns.setVisible(false);
                        NfInput.setEditable(false);
                        NfInput.setVisible(false);
                        Nf.setVisible(false);
                        break;
                    case 5:
                        NInput.setVisible(true);
                        N.setVisible(true);
                        NInput.setEditable(true);
                        NsInput.setEditable(false);
                        NsInput.setVisible(false);
                        Ns.setVisible(false);
                        NfInput.setEditable(false);
                        NfInput.setVisible(false);
                        Nf.setVisible(false);
                        break;
                    case 6:
                        NInput.setVisible(false);
                        N.setVisible(false);
                        NInput.setEditable(false);
                        NsInput.setEditable(true);
                        NsInput.setVisible(true);
                        Ns.setVisible(true);
                        NfInput.setEditable(true);
                        NfInput.setVisible(true);
                        Nf.setVisible(true);
                        break;
                    case 7:
                        NInput.setVisible(true);
                        N.setVisible(true);
                        NInput.setEditable(true);
                        NsInput.setEditable(false);
                        NsInput.setVisible(false);
                        Ns.setVisible(false);
                        NfInput.setEditable(false);
                        NfInput.setVisible(false);
                        Nf.setVisible(false);
                        break;
                    case 8:
                        NInput.setVisible(true);
                        N.setVisible(true);
                        NInput.setEditable(true);
                        NsInput.setEditable(false);
                        NsInput.setVisible(false);
                        Ns.setVisible(false);
                        NfInput.setEditable(false);
                        NfInput.setVisible(false);
                        Nf.setVisible(false);
                        break;
                    case 9:
                        NInput.setVisible(false);
                        N.setVisible(false);
                        NInput.setEditable(false);
                        NsInput.setEditable(true);
                        NsInput.setVisible(true);
                        Ns.setVisible(true);
                        NfInput.setEditable(true);
                        NfInput.setVisible(true);
                        Nf.setVisible(true);
                        break;
                }


            }
        });
        pane.prefHeightProperty().bind(gridPane.heightProperty());
        pane.prefWidthProperty().bind(gridPane.widthProperty());
    }

    public void RungeKutta() {
        yAxis.setForceZeroInRange(false);
        double x0 = Double.parseDouble(x0Input.getText()), X = Double.parseDouble(XInput.getText());
        double y0 = Double.parseDouble(y0Input.getText());
        double c = getConst(x0, y0);
        int N = Integer.parseInt(NInput.getText());
        double[] x = new double[N], y = new double[N];
        x[0] = x0;
        y[0] = y0;
        int index = -1;
        double k1, k2, k3, k4;
        double disc = 2 - Math.log(Math.E - 1);
        double h = (X - x0)/N;
        for (int i = 1; i < N; i++) {
            x[i] = x[i-1] + h;
        }
        for (int i = 1; i < N; i++){
            if (x[i-1] >= disc - h/1.1 && x[i-1] <= disc + h/1.1) {
                index = i;
                i++;
                y[i] = getExact(x[i], c);
            } else {
                k1 = h*(Math.pow(y[i - 1], 2) * Math.pow(Math.E, x[i - 1]) - 2 * y[i - 1]);
                k2 = h*(Math.pow(y[i - 1] + k1/2, 2) * Math.pow(Math.E, x[i - 1] + h/2) - 2 * (y[i - 1] + k1/2));
                k3 = h*(Math.pow(y[i - 1] + k2/2, 2) * Math.pow(Math.E, x[i - 1] + h/2) - 2 * (y[i - 1] + k2/2));
                k4 = h*(Math.pow(y[i - 1] + k3, 2) * Math.pow(Math.E, x[i - 1] + h) - 2 * (y[i - 1] + k3));
                y[i] = y[i - 1] + (k1 + 2*k2 + 2*k3 + k4) / 6;
            }
        }
        XYChart.Series<Number, Number> series = new XYChart.Series<>();
        series.setName("Runge-Kutta");
        for (int i = 0; i < x.length; i++){
            if (i != index) {
                series.getData().add(new XYChart.Data<>(x[i], y[i]));
            }
        }
        chart.getData().add(series);
    }

    public void EulerImproved() {
        yAxis.setForceZeroInRange(false);
        double x0 = Double.parseDouble(x0Input.getText()), X = Double.parseDouble(XInput.getText());
        double y0 = Double.parseDouble(y0Input.getText());
        double c = getConst(x0, y0);
        int N = Integer.parseInt(NInput.getText());
        double[] x = new double[N], y = new double[N];
        x[0] = x0;
        y[0] = y0;
        int index = -1;
        double f = 0;
        double disc = 2 - Math.log(Math.E - 1);
        double h = (X - x0)/N;
        for (int i = 1; i < N; i++) {
            x[i] = x[i-1] + h;
        }
        for (int i = 1; i < N; i++){
            if (x[i-1] >= disc - h/1.1 && x[i-1] <= disc + h/1.1) {
                index = i;
                i++;
                y[i] = getExact(x[i], c);
            } else {
                f = (Math.pow(y[i - 1], 2) * Math.pow(Math.E, x[i - 1]) - 2 * y[i - 1]);
                y[i] = y[i - 1] + h * (Math.pow(y[i - 1] + h/2*f, 2) * Math.pow(Math.E, x[i - 1] + h/2) - 2 * (y[i - 1] + h/2*f));
            }
        }
        XYChart.Series<Number, Number> series = new XYChart.Series<>();
        series.setName("Euler improved");
        for (int i = 0; i < x.length; i++){
            if (i != index) {
                series.getData().add(new XYChart.Data<>(x[i], y[i]));
            }
        }
        chart.getData().add(series);
    }


    public void Euler() {
        yAxis.setForceZeroInRange(false);
        double x0 = Double.parseDouble(x0Input.getText()), X = Double.parseDouble(XInput.getText());
        double y0 = Double.parseDouble(y0Input.getText());
        double c = getConst(x0, y0);
        int N = Integer.parseInt(NInput.getText());
        double[] x = new double[N], y = new double[N];
        x[0] = x0;
        y[0] = y0;
        int index = -1;
        double disc = 2 - Math.log(Math.E - 1);
        double h = (X - x0)/N;
        for (int i = 1; i < N; i++) {
            x[i] = x[i-1] + h;
        }
        for (int i = 1; i < N; i++){
            if (x[i-1] >= disc - h/1.1 && x[i-1] <= disc + h/1.1) {
                index = i;
                i++;
                y[i] = getExact(x[i], c);
           } else {
                y[i] = y[i - 1] + h * (Math.pow(y[i - 1], 2) * Math.pow(Math.E, x[i - 1]) - 2 * y[i - 1]);
            }
        }
        XYChart.Series<Number, Number> series = new XYChart.Series<>();
        series.setName("Euler");
        for (int i = 0; i < x.length; i++){
            if (i != index) {
                series.getData().add(new XYChart.Data<>(x[i], y[i]));
            }
        }
        chart.getData().add(series);
    }

    public void EulerLocalError() {
        double x0 = Double.parseDouble(x0Input.getText()), X = Double.parseDouble(XInput.getText());
        double y0 = Double.parseDouble(y0Input.getText());
        double c = getConst(x0, y0);
        int N = Integer.parseInt(NInput.getText());
        double[] x = new double[N], y = new double[N], yExact = new double[N], error = new double[N];
        x[0] = x0;
        y[0] = y0;
        int index = -1;
        double disc = 2 - Math.log(Math.E - 1);
        double h = (X - x0)/N;
        for (int i = 1; i < N; i++) {
            x[i] = x[i-1] + h;
        }
        for (int i = 1; i < N; i++){
            if (x[i-1] >= disc - h/1.1 && x[i-1] <= disc + h/1.1) {
                index = i;
                i++;
                y[i] = getExact(x[i], c);
            } else {
                y[i] = y[i - 1] + h * (Math.pow(y[i - 1], 2) * Math.pow(Math.E, x[i - 1]) - 2 * y[i - 1]);
            }
        }
        for (int i = 0; i < N; i++){
            yExact[i] = Math.pow(Math.E, 2 - x[i]) / (Math.pow(Math.E, x[i]) - Math.pow(Math.E, x[i] + 1) + Math.pow(Math.E, 2));
            if (i != index - 1) {
                error[i] = Math.abs(yExact[i] - y[i]);
            }
        }

        XYChart.Series<Number, Number> series = new XYChart.Series<>();
        series.setName("Euler local error");
        for (int i = 0; i < x.length; i++){
            series.getData().add(new XYChart.Data<>(x[i], error[i]));
        }
        chart.getData().add(series);

    }

    public void Exact() {
        double x0 = Double.parseDouble(x0Input.getText()), X = Double.parseDouble(XInput.getText());
        double y0 = Double.parseDouble(y0Input.getText());
        int N = Integer.parseInt(NInput.getText());
        double[] x = new double[N], y = new double[N];
        x[0] = x0;
        double c = getConst(x0, y0);
        double disc = 2 - Math.log(Math.E - 1);
        int index = -1;
        double h = (X - x0)/N;
        for (int i = 1; i < N; i++) {
            x[i] = x[i-1] + h;
        }
        for (int i = 0; i < N; i++){
            if (x[i] >= disc - h/1.1 && x[i] <= disc + h/1.1) {
                index = i;
                i++;
            }
            y[i] = getExact(x[i], c);
        }
        XYChart.Series<Number, Number> series = new XYChart.Series<>();
        series.setName("Exact");
        for (int i = 0; i < x.length; i++){
            if (i != index) {
                series.getData().add(new XYChart.Data<>(x[i], y[i]));
            }
        }
        chart.getData().add(series);
    }

    public void EulerTotalError() {
        double x0 = Double.parseDouble(x0Input.getText()), X = Double.parseDouble(XInput.getText());
        double y0 = Double.parseDouble(y0Input.getText());
        double c = getConst(x0, y0);
        int Ns = Integer.parseInt(NsInput.getText()), Nf = Integer.parseInt(NfInput.getText());
        double[] x = new double[Nf], y = new double[Nf], yExact = new double[Nf], error = new double[Nf];
        double h;
        double max;
        int index = -1;
        x[0] = x0;
        y[0] = y0;
        double disc = 2 - Math.log(Math.E - 1);
        for (int i = Ns; i <= Nf; i++) {
            h = (X - x0) / i;
            max = 0;
            yExact[0] = y[0];
            for (int k = 1; k < i; k++) {
                x[k] = x[k-1] + h;
                if (x[k] >= disc - h/1.1 && x[k] <= disc + h/1.1) {
                    index = k;
                }
                yExact[k] = getExact(x[k], c);
            }
            yExact[index] = 0;
            for (int k = 1; k < i; k++){
                if (x[k-1] >= disc - h/1.1 && x[k-1] <= disc + h/1.1) {
                    k++;
                    y[k] = getExact(x[k], c);
                } else {
                    y[k] = y[k - 1] + h * (Math.pow(y[k - 1], 2) * Math.pow(Math.E, x[k - 1]) - 2 * y[k - 1]);
                    if (max < Math.abs(yExact[k]-y[k])) { max = Math.abs(yExact[k]-y[k]); }
                }
            }
            error[i-Ns] = max;
        }
        XYChart.Series<Number, Number> series = new XYChart.Series<>();
        series.setName("Euler total error");
        for (int i = Ns; i <= Nf; i++){
            series.getData().add(new XYChart.Data<>(i, error[i-Ns]));
        }
        chart.getData().add(series);

    }

    public void EulerImprovedLocalError() {
        double x0 = Double.parseDouble(x0Input.getText()), X = Double.parseDouble(XInput.getText());
        double y0 = Double.parseDouble(y0Input.getText());
        double c = getConst(x0, y0);
        int N = Integer.parseInt(NInput.getText());
        double[] x = new double[N], y = new double[N], yExact = new double[N], error = new double[N];
        x[0] = x0;
        y[0] = y0;
        int index = -1;
        double disc = 2 - Math.log(Math.E - 1);
        double h = (X - x0)/N;
        double f;
        for (int i = 1; i < N; i++) {
            x[i] = x[i-1] + h;
        }
        for (int i = 1; i < N; i++){
            if (x[i-1] >= disc - h/1.1 && x[i-1] <= disc + h/1.1) {
                index = i;
                i++;
                y[i] = getExact(x[i], c);
            } else {
                f = (Math.pow(y[i - 1], 2) * Math.pow(Math.E, x[i - 1]) - 2 * y[i - 1]);
                y[i] = y[i - 1] + h * (Math.pow(y[i - 1] + h/2*f, 2) * Math.pow(Math.E, x[i - 1] + h/2) - 2 * (y[i - 1] + h/2*f));
            }
        }
        for (int i = 0; i < N; i++){
            yExact[i] = getExact(x[i], c);
            if (i != index - 1) {
                error[i] = Math.abs(yExact[i] - y[i]);
            }
        }
        XYChart.Series<Number, Number> series = new XYChart.Series<>();
        series.setName("Euler improved local error");
        for (int i = 0; i < x.length; i++){
            series.getData().add(new XYChart.Data<>(x[i], error[i]));
        }
        chart.getData().add(series);

    }

    public void EulerImprovedTotalError() {
        double x0 = Double.parseDouble(x0Input.getText()), X = Double.parseDouble(XInput.getText());
        double y0 = Double.parseDouble(y0Input.getText());
        double c = getConst(x0, y0);
        int Ns = Integer.parseInt(NsInput.getText()), Nf = Integer.parseInt(NfInput.getText());
        double[] x = new double[Nf], y = new double[Nf], yExact = new double[Nf], error = new double[Nf];
        double h;
        double f;
        int index = -1;
        double max;
        x[0] = x0;
        y[0] = y0;
        double disc = 2 - Math.log(Math.E - 1);
        for (int i = Ns; i <= Nf; i++) {
            h = (X - x0) / i;
            max = 0;
            yExact[0] = y[0];
            for (int k = 1; k < i; k++) {
                x[k] = x[k-1] + h;
                if (x[k] >= disc - h/1.1 && x[k] <= disc + h/1.1) {
                    index = k;
                }
                yExact[k] = getExact(x[k], c);
            }
            yExact[index] = 0;
            for (int k = 1; k < i; k++){
                if (x[k-1] >= disc - h/1.1 && x[k-1] <= disc + h/1.1) {
                    k++;
                    y[k] = getExact(x[k], c);
                } else {
                    f = (Math.pow(y[k - 1], 2) * Math.pow(Math.E, x[k - 1]) - 2 * y[k - 1]);
                    y[k] = y[k - 1] + h * (Math.pow(y[k - 1] + h/2*f, 2) * Math.pow(Math.E, x[k - 1] + h/2) - 2 * (y[k - 1] + h/2*f));
                    if (max < Math.abs(yExact[k]-y[k])) { max = Math.abs(yExact[k]-y[k]); }
                }
            }
            error[i-Ns] = max;
        }
        XYChart.Series<Number, Number> series = new XYChart.Series<>();
        series.setName("Euler improved total error");
        for (int i = Ns; i <= Nf; i++){
            series.getData().add(new XYChart.Data<>(i, error[i-Ns]));
        }
        chart.getData().add(series);

    }

    public void RungeKuttaLocalError() {
        double x0 = Double.parseDouble(x0Input.getText()), X = Double.parseDouble(XInput.getText());
        double y0 = Double.parseDouble(y0Input.getText());
        double c = getConst(x0, y0);
        int N = Integer.parseInt(NInput.getText());
        double[] x = new double[N], y = new double[N], yExact = new double[N], error = new double[N];
        x[0] = x0;
        y[0] = y0;
        int index = -1;
        double disc = 2 - Math.log(Math.E - 1);
        double h = (X - x0)/N;
        double k1, k2, k3, k4;
        for (int i = 1; i < N; i++) {
            x[i] = x[i-1] + h;
        }
        for (int i = 1; i < N; i++){
            if (x[i-1] >= disc - h/1.1 && x[i-1] <= disc + h/1.1) {
                index = i;
                i++;
                y[i] = getExact(x[i], c);
            } else {
                k1 = h*(Math.pow(y[i - 1], 2) * Math.pow(Math.E, x[i - 1]) - 2 * y[i - 1]);
                k2 = h*(Math.pow(y[i - 1] + k1/2, 2) * Math.pow(Math.E, x[i - 1] + h/2) - 2 * (y[i - 1] + k1/2));
                k3 = h*(Math.pow(y[i - 1] + k2/2, 2) * Math.pow(Math.E, x[i - 1] + h/2) - 2 * (y[i - 1] + k2/2));
                k4 = h*(Math.pow(y[i - 1] + k3, 2) * Math.pow(Math.E, x[i - 1] + h) - 2 * (y[i - 1] + k3));
                y[i] = y[i - 1] + (k1 + 2*k2 + 2*k3 + k4) / 6;
            }
        }
        for (int i = 0; i < N; i++){
            yExact[i] = getExact(x[i], c);
            if (i != index - 1) {
                error[i] = Math.abs(yExact[i] - y[i]);
            }
        }
        XYChart.Series<Number, Number> series = new XYChart.Series<>();
        series.setName("Runge-Kutta local error");
        for (int i = 0; i < x.length; i++){
            series.getData().add(new XYChart.Data<>(x[i], error[i]));
        }
        chart.getData().add(series);

    }

    public void RungeKuttaTotalError() {
        double x0 = Double.parseDouble(x0Input.getText()), X = Double.parseDouble(XInput.getText());
        double y0 = Double.parseDouble(y0Input.getText());
        double c = getConst(x0, y0);
        int Ns = Integer.parseInt(NsInput.getText()), Nf = Integer.parseInt(NfInput.getText());
        double[] x = new double[Nf], y = new double[Nf], yExact = new double[Nf], error = new double[Nf];
        double h;
        double k1, k2, k3, k4;
        double max;
        int index = -1;
        x[0] = x0;
        y[0] = y0;
        double disc = 2 - Math.log(Math.E - 1);
        for (int i = Ns; i <= Nf; i++) {
            h = (X - x0) / i;
            max = 0;
            yExact[0] = y[0];
            for (int k = 1; k < i; k++) {
                x[k] = x[k-1] + h;
                if (x[k] >= disc - h/1.1 && x[k] <= disc + h/1.1) {
                    index = k;
                }
                    yExact[k] = getExact(x[k], c);
            }
            yExact[index] = 0;
            for (int k = 1; k < i; k++){
                if (x[k-1] >= disc - h/1.1 && x[k-1] <= disc + h/1.1) {
                    k++;
                    y[k] = getExact(x[k], c);
                } else {
                    k1 = h*(Math.pow(y[k - 1], 2) * Math.pow(Math.E, x[k - 1]) - 2 * y[k - 1]);
                    k2 = h*(Math.pow(y[k - 1] + k1/2, 2) * Math.pow(Math.E, x[k - 1] + h/2) - 2 * (y[k - 1] + k1/2));
                    k3 = h*(Math.pow(y[k - 1] + k2/2, 2) * Math.pow(Math.E, x[k - 1] + h/2) - 2 * (y[k - 1] + k2/2));
                    k4 = h*(Math.pow(y[k - 1] + k3, 2) * Math.pow(Math.E, x[k - 1] + h) - 2 * (y[k - 1] + k3));
                    y[k] = y[k - 1] + (k1 + 2*k2 + 2*k3 + k4) / 6;
                    if (max < Math.abs(yExact[k]-y[k])) { max = Math.abs(yExact[k]-y[k]); }
                }
            }
            error[i-Ns] = max;
        }
        XYChart.Series<Number, Number> series = new XYChart.Series<>();
        series.setName("Runge-Kutta total error");
        for (int i = Ns; i <= Nf; i++){
            series.getData().add(new XYChart.Data<>(i, error[i-Ns]));
        }
        chart.getData().add(series);

    }

    public double getConst(double x, double y) {
        return (1 - y * Math.pow(Math.E, x)) / (y * Math.pow(Math.E, 2*x));
    }

    public double getExact(double x, double c) {
        return 1 / (Math.pow(Math.E, x) * (c * Math.pow(Math.E, x) + 1));
    }



    public void plot() {
        Alert alert = new Alert(Alert.AlertType.ERROR);
        alert.setTitle("ERROR");

        alert.setContentText("You have typed wrong input");
        if (!(x0Input.getText().matches("-?\\d+(\\.\\d+)?") && XInput.getText().matches("-?\\d+(\\.\\d+)?") && y0Input.getText().matches("-?\\d+(\\.\\d+)?"))) {
            alert.showAndWait();
        } else {
            switch (choiceBox.getValue()) {
                case "Exact":
                    if (!NInput.getText().matches("\\d+")){
                        alert.showAndWait();
                        break;
                    }
                    Exact();
                    break;
                case "Euler":
                    if (!NInput.getText().matches("\\d+")){
                        alert.showAndWait();
                        break;
                    }
                    Euler();
                    break;
                case "Euler local error":
                    if (!NInput.getText().matches("\\d+")){
                        alert.showAndWait();
                        break;
                    }
                    EulerLocalError();
                    break;
                case "Euler total error":
                    if (!(NsInput.getText().matches("\\d+") && NfInput.getText().matches("\\d+") && (Integer.parseInt(NsInput.getText()) < Integer.parseInt(NfInput.getText())))){
                        alert.showAndWait();
                        break;
                    }
                    EulerTotalError();
                    break;
                case "Euler improved":
                    if (!NInput.getText().matches("\\d+")){
                        alert.showAndWait();
                        break;
                    }
                    EulerImproved();
                    break;
                case "Euler improved local error":
                    if (!NInput.getText().matches("\\d+")){
                        alert.showAndWait();
                        break;
                    }
                    EulerImprovedLocalError();
                    break;
                case "Euler improved total error":
                    if (!(NsInput.getText().matches("\\d+") && NfInput.getText().matches("\\d+") && (Integer.parseInt(NsInput.getText()) < Integer.parseInt(NfInput.getText())))){
                        alert.showAndWait();
                        break;
                    }
                    EulerImprovedTotalError();
                    break;
                case "Runge-Kutta":
                    if (!NInput.getText().matches("\\d+")){
                        alert.showAndWait();
                        break;
                    }
                    RungeKutta();
                    break;
                case "Runge-Kutta local error":
                    if (!NInput.getText().matches("\\d+")){
                        alert.showAndWait();
                        break;
                    }
                    RungeKuttaLocalError();
                    break;
                case "Runge-Kutta total error":
                    if (!(NsInput.getText().matches("\\d+") && NfInput.getText().matches("\\d+") && (Integer.parseInt(NsInput.getText()) < Integer.parseInt(NfInput.getText())))){
                        alert.showAndWait();
                        break;
                    }
                    RungeKuttaTotalError();
                    break;
            }
        }

    }

    public void clear() {
        chart.getData().clear();
    }






    public void quit() {
        Stage stage = (Stage) quitButton.getScene().getWindow();
        stage.close();
    }

}
