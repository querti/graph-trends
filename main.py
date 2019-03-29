
import csv
import sys
import matplotlib.pyplot as plt
import numpy as np


class GraphPlotter:

    def __init__(self, paths, datas, ranges, lengths, show_trend = False):

        self.paths = paths
        self.datas = datas
        self.ranges = ranges
        self.lengths = lengths
        self.show_trend = show_trend

    def make_graph(self):

        plots = []
        it = 0
        trend_equations = []
        for path in self.paths:

            start_x, end_x = self.ranges[it]
            start_x = float(start_x)
            end_x = float(end_x)
            

            position_list = []
            load_list = []
            deformation_list = []
            deformation_corrected_list = []
            stress_list = []

            suffix = path.split('/')[-1][4]
            height = float(self.datas[suffix][0])
            surface_area = float(self.datas[suffix][1])
            name = path.split('/')[-1].split('.')[0]

            with open(path) as csv_f:
                csv_contracts = csv.DictReader(csv_f, delimiter=',')
                for row in csv_contracts:
                    position_list.append(float(row['Position(Tah-Tlak:Position) (mm)']))
                    load_list.append(float(row['Load(Tah-Tlak:Load) (kN)']))
                    deformation = (float(row['Position(Tah-Tlak:Position) (mm)']) * -1) / height
                    deformation_corrected = deformation - 1.7787
                    stress = (float(row['Load(Tah-Tlak:Load) (kN)']) * (-1000)) / surface_area

                    deformation_list.append(deformation)
                    deformation_corrected_list.append(deformation_corrected)
                    stress_list.append(stress)

            plt.plot(deformation_corrected_list, stress_list, label=eval("path.split('/')[-1]"))

            if self.show_trend:
                self.find_best_regression(start_x, end_x, deformation_corrected_list, stress_list, name)
            else:
                with open('results/DENSE CONTROL/' + name+'-regressions.csv', 'r') as f:
                    csv_contracts = csv.DictReader(f, delimiter=',')
                    for row in csv_contracts:
                        if int(row['regression_length']) == self.lengths[it]:
                            start_x_p = float(row['start_x'])
                            end_x_p = float(row['end_x'])
                            a = float(row['a'])
                            b = float(row['b'])

                x = np.array([start_x_p, end_x_p])  
                y = eval('a + b*x')
                trend_equations.append(name+ ': y = {}*x + {}'.format(b,a))
                plt.plot(x, y) 
            it += 1
        
        with open('results/DENSE CONTROL/trends.txt', 'w+') as f:
            for i in trend_equations:
                f.write(i+'\n')
        plt.legend()
        plt.title('DENSE CONTROL')
        plt.show()

    def find_best_regression(self, start_x, end_x, x, y, name):
        part_x = []
        part_y = []
        best_regressions = []

        for i in range(len(x)):
            if x[i] >= start_x and x[i] <= end_x:
                part_x.append(x[i])
                part_y.append(y[i])
        
        for i in range(len(part_x)):
            if i <= 2:
                continue
            best_a, best_b, best_err, best_x1, best_x2 = self.calculate_len_regressions(part_x, part_y, i)
            best_regressions.append((i, best_a, best_b, best_err, best_x1, best_x2))
        
        best = sorted(best_regressions, key=lambda x: x[3])

        with open('results/DENSE CONTROL/' + name+'-regressions.csv', 'w+') as f:
            f.write('regression_length,a,b,error,start_x,end_x\n')
            for i in best:
                f.write('{},{},{},{},{},{}\n'.format(i[0],i[1],i[2],i[3],i[4],i[5]))
    
    def calculate_len_regressions(self, x, y, length):
        best_a = 0
        best_b = 0
        best_err = 100000000000000
        best_x1 = 0
        best_x2 = 0
        fields_x = []
        fields_y = []

        for i in range(len(x) - length):
            field_x = []
            field_y = []
            for j in range(length):
                field_x.append(x[i+j])
                field_y.append(y[i+j])
            
            fields_x.append(field_x)
            fields_y.append(field_y)
        
        for i in range(len(fields_x)):
            a, b = self.calculate_partial_regression(0, 0, fields_x[i], fields_y[i], True)
            y_estimates = []
            for j in fields_x[i]:
                y_estimates.append(a + b*j)
            error = self.calculate_mean_square_error(fields_y[i], y_estimates)
            if error < best_err:
                best_a = a
                best_b = b
                best_err = error
                best_x1 = fields_x[i][0]
                best_x2 = fields_x[i][-1]
            
        return(best_a, best_b, best_err, best_x1, best_x2)
    
    def calculate_partial_regression(self, start_x, end_x, x, y, whole_plot= False):

        sum_x = 0
        sum_y = 0
        sum_x_sq = 0
        sum_xy = 0
        counter = 0
        for i in range(0, len(x)):
            if (x[i] > start_x and x[i] < end_x) or whole_plot:
                sum_x += x[i]
                sum_y += y[i]
                sum_x_sq += x[i]*x[i]
                sum_xy += x[i]*y[i]
                counter += 1
        
        n = counter
        a = (sum_y*sum_x_sq - sum_x*sum_xy) / (n*sum_x_sq - sum_x*sum_x)
        b = (n*sum_xy - sum_x*sum_y) / (n*sum_x_sq - sum_x*sum_x)

        return (a, b)
    
    def calculate_mean_square_error(self, true_y, est_y):
        summ = 0
        for i in range(len(true_y)):
            summ += (true_y[i] - est_y[i])*(true_y[i] - est_y[i])
        
        return summ / len(true_y)
    
def main():
    """Main function."""
    datas = {}
    with open('datas.csv','r') as f:
        while True:
            line = f.readline()
            if 'height' in line:
                continue
            if line == '':
                break
            i, h, a = line[:-1].split(',')
            datas[i] = (h, a)

    items_dense_control = ['190327_CDHA_cylinders/DENSE CONTROL/Test6/Test6.Stop.csv',
                           '190327_CDHA_cylinders/DENSE CONTROL/Test7/Test7.Stop.csv',
                           '190327_CDHA_cylinders/DENSE CONTROL/Test8/Test8.Stop.csv',
                           '190327_CDHA_cylinders/DENSE CONTROL/Test9/Test9.Stop.csv']
    
    lengths_dense_control = [150,
                             150,
                             200,
                             250]
    ranges_dense_control = [(-0.363,-0.3559),
                            (-0.378,-0.371),
                            (-0.372,-0.359),
                            (-0.165,-0.157)]
    
    items_dense_pcl1 = ['190327_CDHA_cylinders/DENSE PCL 1/Test15/Test15.Stop.csv',
                        '190327_CDHA_cylinders/DENSE PCL 1/Test16/Test16.Stop.csv',
                        '190327_CDHA_cylinders/DENSE PCL 1/Test17/Test17.Stop.csv',
                        '190327_CDHA_cylinders/DENSE PCL 1/Test18/Test18.Stop.csv',
                        '190327_CDHA_cylinders/DENSE PCL 1/Test19/Test19.Stop.csv']
    
    lines_dense_pcl1 = [(-0.106,-0.103),
                        (-0.127,-0.122),
                        (-0.040,-0.034),
                        (-0.125,-0.119),
                        (-0.133,-0.127)]

    items_dense_plu1 = ['190327_CDHA_cylinders/DENSE PLU 1/Test10/Test10.Stop.csv',
                        '190327_CDHA_cylinders/DENSE PLU 1/Test11/Test11.Stop.csv',
                        '190327_CDHA_cylinders/DENSE PLU 1/Test12/Test12.Stop.csv',
                        '190327_CDHA_cylinders/DENSE PLU 1/Test13/Test13.Stop.csv',
                        '190327_CDHA_cylinders/DENSE PLU 1/Test14/Test14.Stop.csv']
    lines_dense_plu1 = [(-0.060,-0.056),
                        (-0.036,-0.030),
                        (-0.077,-0.069),
                        (-0.009,-0.003),
                        (0.043,0.048)]

    plotter = GraphPlotter(items_dense_control, datas, ranges_dense_control, lengths_dense_control, False)
    plotter.make_graph()

if __name__ == "__main__":
    # cProfile.run('main()')
    main()