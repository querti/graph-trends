
import csv
import sys
import matplotlib.pyplot as plt
import numpy as np


class GraphPlotter:

    def __init__(self, paths, datas, ranges, lengths, show_trend = False, skip = False):

        self.paths = paths
        self.datas = datas
        self.ranges = ranges
        self.lengths = lengths
        self.show_trend = show_trend
        self.skip = skip

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

            if self.show_trend and not self.skip:
                self.find_best_regression(start_x, end_x, deformation_corrected_list, stress_list, name)
            elif not self.show_trend and not self.skip:
                with open('results/POROUS PCL 1/' + name+'-regressions.csv', 'r') as f:
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
        
        with open('results/POROUS PCL 1/trends.txt', 'w+') as f:
            for i in trend_equations:
                f.write(i+'\n')
        if self.skip or not self.show_trend:
            plt.legend()
        plt.title('POROUS PCL 1')
        plt.xlabel('Deformation [mm/mm]')
        plt.ylabel('Stress [MPa]')
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

        with open('results/POROUS PCL 1/' + name+'-regressions.csv', 'w+') as f:
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
    
    lengths_dense_control = [100,
                             150,
                             250,
                             200]
    ranges_dense_control = [(-0.363,-0.3559),
                            (-0.378,-0.371),
                            (-0.372,-0.359),
                            (-0.165,-0.157)]
    
    items_dense_pcl1 = ['190327_CDHA_cylinders/DENSE PCL 1/Test15/Test15.Stop.csv',
                        '190327_CDHA_cylinders/DENSE PCL 1/Test16/Test16.Stop.csv',
                        '190327_CDHA_cylinders/DENSE PCL 1/Test17/Test17.Stop.csv',
                        '190327_CDHA_cylinders/DENSE PCL 1/Test18/Test18.Stop.csv',
                        '190327_CDHA_cylinders/DENSE PCL 1/Test19/Test19.Stop.csv']
    
    lengths_dense_pcl1 = [130,
                          150,
                          250,
                          150,
                          180]
    ranges_dense_pcl1 = [(-0.108,-0.099),
                        (-0.127,-0.101),
                        (-0.043,-0.029),
                        (-0.126,-0.118),
                        (-0.135,-0.121)]

    items_dense_plu1 = ['190327_CDHA_cylinders/DENSE PLU 1/Test10/Test10.Stop.csv',
                        '190327_CDHA_cylinders/DENSE PLU 1/Test11/Test11.Stop.csv',
                        '190327_CDHA_cylinders/DENSE PLU 1/Test12/Test12.Stop.csv',
                        '190327_CDHA_cylinders/DENSE PLU 1/Test13/Test13.Stop.csv',
                        '190327_CDHA_cylinders/DENSE PLU 1/Test14/Test14.Stop.csv']
    lengths_dense_plu1 = [150,
                          150,
                          200,
                          150,
                          150]
    ranges_dense_plu1 = [(-0.070,-0.043),
                        (-0.039,-0.029),
                        (-0.082,-0.069),
                        (-0.020,-0.0001),
                        (0.039,0.051)]

    items_porous_control = ['190327_CDHA_cylinders/POROUS CONTROL/Test1/Test1.Stop.csv',
                            '190327_CDHA_cylinders/POROUS CONTROL/Test2/Test2.Stop.csv',
                            '190327_CDHA_cylinders/POROUS CONTROL/Test3/Test3.Stop.csv',
                            '190327_CDHA_cylinders/POROUS CONTROL/Test4/Test4.Stop.csv',
                            '190327_CDHA_cylinders/POROUS CONTROL/Test5/Test5.Stop.csv']
    lengths_porous_control = [160,
                              150,
                              150,
                              180,
                              150]
    ranges_porous_control = [(0.0015,0.0118),
                             (-0.2032,-0.1973),
                             (-0.0673,-0.0605),
                             (0.0887,0.104),
                             (-0.0346,-0.0247)]
    
    items_porous_pcl1 = ['190327_CDHA_cylinders/POROUS PCL 1/Test45/Test45.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PCL 1/Test46/Test46.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PCL 1/Test47/Test47.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PCL 1/Test48/Test48.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PCL 1/Test49/Test49.Stop.csv']
    lengths_porous_pcl1 = [140,
                           140,
                           200,
                           190,
                           160]
    ranges_porous_pcl1 = [(0.0034,0.0115),
                          (0.051,0.0587),
                          (-0.0227,-0.0106),
                          (0.0633,0.0720),
                          (0.0466,0.0597)]

    items_porous_pcl3 = ['190327_CDHA_cylinders/POROUS PCL 3/Test40/Test40.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PCL 3/Test41/Test41.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PCL 3/Test42/Test42.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PCL 3/Test43/Test43.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PCL 3/Test44/Test44.Stop.csv']
    lengths_porous_pcl3 = [140,
                           160,
                           200,
                           170,
                           160]
    ranges_porous_pcl3 = [(0.0656,0.0746),
                          (0.0463,0.0607),
                          (0.0022,0.0117),
                          (-0.0559,-0.0478),
                          (0.0977,0.1099)]
    
    items_porous_pcl5 = ['190327_CDHA_cylinders/POROUS PCL 5/Test35/Test35.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PCL 5/Test36/Test36.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PCL 5/Test37/Test37.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PCL 5/Test38/Test38.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PCL 5/Test39/Test39.Stop.csv']
    lengths_porous_pcl5 = [140,
                           200,
                           180,
                           170,
                           160]
    ranges_porous_pcl5 = [(-0.1352,-0.1212),
                          (-0.0835,-0.0720),
                          (-0.0956,-0.0846),
                          (-0.0537,-0.0423),
                          (-0.1254,-0.1154)]

    items_porous_plu1 = ['190327_CDHA_cylinders/POROUS PLU 1/Test20/Test20.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PLU 1/Test21/Test21.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PLU 1/Test22/Test22.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PLU 1/Test23/Test23.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PLU 1/Test24/Test24.Stop.csv']
    lengths_porous_plu1 = [140,
                           200,
                           180,
                           170,
                           160]
    ranges_porous_plu1 = [(-0.1044,-0.0982),
                          (-0.1389,-0.1311),
                          (-0.1356,-0.1148),
                          (-0.1765,-0.1677),
                          (-0.1600,-0.1479)]
    
    items_porous_plu3 = ['190327_CDHA_cylinders/POROUS PLU 3/Test30/Test30.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PLU 3/Test31/Test31.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PLU 3/Test32/Test32.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PLU 3/Test33/Test33.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PLU 3/Test34/Test34.Stop.csv']
    lengths_porous_plu3 = [140,
                           200,
                           180,
                           170,
                           160]
    ranges_porous_plu3 = [(-0.1120,-0.1021),
                          (-0.1946,-0.1855),
                          (-0.1313,-0.1214),
                          (-0.1619,-0.1508),
                          (-0.1077,-0.1003)]

    items_porous_plu5 = ['190327_CDHA_cylinders/POROUS PLU 5/Test25/Test25.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PLU 5/Test26/Test26.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PLU 5/Test27/Test27.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PLU 5/Test28/Test28.Stop.csv',
                         '190327_CDHA_cylinders/POROUS PLU 5/Test29/Test29.Stop.csv']
    lengths_porous_plu5 = [140,
                           200,
                           180,
                           170,
                           160]
    ranges_porous_plu5 = [(-0.2387,-0.2309),
                          (-0.1953,-0.1858),
                          (-0.1769,-0.1685),
                          (-0.2126,-0.2007),
                          (-0.2004,-0.1892)]

    plotter = GraphPlotter(items_porous_pcl1, datas, ranges_porous_pcl1, lengths_porous_pcl1, False, False)
    plotter.make_graph()

if __name__ == "__main__":
    # cProfile.run('main()')
    main()