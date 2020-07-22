from scipy.integrate import quad
from numpy import inf, linspace
from matplotlib import pyplot as plt


class MemberFunct:
    def __init__(self,field_name, memb_name,triangle_func=True,evaluate_funct=None,*args):
        self.triangle_func = triangle_func
        self.field_name = field_name
        self.memb_name = memb_name
        self.triangle_func = triangle_func
        self.args = args
        if evaluate_funct is None:
            if triangle_func:
                a, b, c = args
                self.evaluate = self.triangle(a, b, c)
                self.points = [(a,0),(b,1),(c,0)]
            else:
                a, b, c, d = args
                self.evaluate = self.trapezoid(a, b, c, d)
                if c is inf:
                    self.points = [(a, 0), (b, 1), (c, 1), (d, 1)]
                else:
                    self.points = [(a, 0), (b, 1), (c, 1), (d, 0)]
        else:
            self.evaluate = evaluate_funct

    def _update_points(self, y):
        '''
        Recalcula los puntos de cambio de recta, dado un corte 'y'
        paralelo al eje de las x.
        :param y: valor de la recta del corte a insertar.
        :return:
        '''
        if self.triangle_func:
            point1 = self.points[0]
            point2 = self.points[1]
            result_1 = self._get_x_from_y(point1, point2, y)

            point1 = self.points[1]
            point2 = self.points[2]
            result_2 = self._get_x_from_y(point1, point2, y)
        else:
            point1 = self.points[0]
            point2 = self.points[1]
            result_1 = self._get_x_from_y(point1, point2, y)

            point1 = self.points[2]
            point2 = self.points[3]

            if point1[0] != inf and point2[0] != inf:
                result_2 = self._get_x_from_y(point1, point2, y)
            else:
                result_2 = "Error"
        self.points = [(xi, yi) for xi, yi in self.points if yi < y]
        self.points.append(result_1)
        if result_2 != 'Error':
            self.points.append(result_2)
        else:
            self.points.append((inf, y))
            self.points.append((inf, y))
        self.points.sort()

    def _get_x_from_y(self, point1, point2, y):
        '''
        Dado el valor de la imagen 'y', trata de calcular su valor 'x'
        :param point1: Un punto de la recta
        :param point2: Otro punto sobre la recta
        :param y: El valor al cual se le busca la x.
        :return: la x
        '''
        den = (point1[0] - point2[0])
        if den == 0:
            den = -1e-2
        m = (point1[1] - point2[1]) / den
        n = point1[1] - m * point1[0]
        inverse = lambda y, m, n: (y - n) / m if m != 0 else 0
        result1 = (inverse(y, m, n), y)
        return result1

    def __add__(self, other):
        evaluate = lambda x: max(self.evaluate(x), other.evaluate(x))
        new_memb = MemberFunct(self.field_name,self.memb_name,self.triangle_func, evaluate)
        new_memb.points = self.points + other.points

        return new_memb


    @staticmethod
    def triangle(a, b, c):
        def evaluate(x):
            return max(min((x - a) / (b - a), (c - x) / (c - b)), 0)

        return evaluate

    @staticmethod
    def trapezoid(a, b, c, d):
        def evaluate(x):
            if b-a != 0:
                return max(min((x - a) / (b - a), 1, (d - x) / (d - c)), 0)
            return max(min(1, (d - x) / (d - c)), 0)

        return evaluate

    def plot_aprox_function(self,**args):
        '''
        Se realiza un ploteo discreto de la funcion de membressía
        :param args: puede recibir como parámetros a:
        - start: entero que indica desde donde comenzara plotear.
            Valor por defecto: 0
        - stop: entero que indica hasta donde comenzar a plotear.
            Valor por defecto:el maximo de la lista self.points distinto de inf
        - num: cantidad de numeros a evaluar sore la función.
            Valor por defecto: (stop-start)*1e3
        :return:
        '''
        self.points.sort()
        mina = self.points[0][0]
        maxa = 0

        for i in range(len(self.points)-1, 0, -1):
            maxa = self.points[i][0]
            if maxa is not inf:
                maxa += 10
                break

        start = mina if 'start' not in args else args['start']
        stop = maxa if 'stop' not in args else args['stop']
        num = (maxa-mina)*1e3 if 'num' not in args else args['num']

        pointsX = linspace(start=start, stop=stop, num=num)
        pointsY = [self.evaluate(x) for x in pointsX]

        plt.plot(pointsX, pointsY)  # plot x and y using default line style and color
        plt.show()


class DiffSystem:

    def __init__(self, mapFieldToMemberFunctName: dict):
        '''
        :param mapFieldToMemberFunctName: diccionario que mapea de ('Field','MemberFuctionName') to
         MemberFunction. Ejemlo: ('Temperatura','Caliente') -> <memberFunction>
        '''
        self.rules = []
        # self.member_functions = memberfuncts
        # self.deffts = deffuzzificators
        self.mapFieldToMemberFunctName = mapFieldToMemberFunctName
        self.input_data = None
        self.defuzzy = {
            'mom': self.mom,
            'som': self.som,
            'lom': self.lom,
            'bisec': self.bisec,
            'centroid': self.centroid
        }
        self.agregation_methods = {
            'mamdani': self.mamdani,
            'larsen': self.larsen
        }

    def add_rule(self, if_body, then_body):
        '''
        Añade reglas al sistema
        :param if_body: composiciones de fuzzy_and, fuzzy_or, fuzzy_not.
        :param then_body: str con el siguiente formato: 'field member_function_name'. Ejemplo: 'Configuracion A'
        :return:
        '''
        self.rules.append((if_body, then_body))

    def evaluate_rules(self, input_data):
        self.set_input_data(input_data)
        return [rule() for rule, _ in self.rules]

    def get_output(self, input_data, agregation_fuction: str, defuzzyficator_function: str,plot=True):
        '''
        :param input_data: diccionario que contiene como llave el nombre de cada campo, y como valor, el valor de ese
         campo obtenido en la medición.

        :param agregation_fuction: str con el nombre del metodo de agregacion q desea utilizar. Disponibles:
        mamdani y larsen.
        :param defuzzyficator_function: nombre del defusificador a utilizar. Disponibles: centroid, mom, som, bisec, lom.
        :return: resultado de la inferencia y la lista de las evaluaciones de las reglas
        '''
        self.set_input_data(input_data)
        agregation_fuction = self.agregation_methods[agregation_fuction]
        defuzzyficator_function = self.defuzzy[defuzzyficator_function]

        rules_evaluations = self.get_rules_evaluations()
        final_member_f = None
        for if_evaluation, then_body in rules_evaluations:

            field, member_f_name = then_body.split()
            member_f = self.get_member_func(field, member_f_name)
            member_f = agregation_fuction(if_evaluation, member_f)

            if final_member_f is None:
                final_member_f = member_f
            else:
                final_member_f += member_f

        result = defuzzyficator_function(final_member_f)
        if plot:
            final_member_f.plot_aprox_function()
        return result, rules_evaluations

    def get_rules_evaluations(self):
        rule_dict = {then_body: [] for if_body,then_body in self.rules}
        for if_body, then_body in self.rules:
            rule_dict[then_body].append(if_body())
        result = []
        for key in rule_dict:
            result.append((max(rule_dict[key]), key))

        return result

    def get_member_func(self,field_name,member_f_name):
        '''
        Devuelve una instancia de la función de membresía member_f_name
        :param field_name: nombre del parámetro que describe la función de membresía
        :param member_f_name: nombre de la función de membresía
        :return: Una nueva instancia de la función de membresía.
        '''
        member_f = self.mapFieldToMemberFunctName[(field_name,member_f_name)]
        new_member_f = MemberFunct(field_name,member_f_name,member_f.triangle_func,None,*member_f.args)
        return new_member_f

    def require_call_set_input_method(self):
        if self.input_data is None:
            print("First, you need to call set_input_data method")
            return

    def fuzzy_and(self,param1,param2):
        fetch_param_result = self._fetch_param_result

        def inside():
            result_param1 = fetch_param_result(param1)
            result_param2 = fetch_param_result(param2)

            return min(result_param1, result_param2)

        return inside

    def fuzzy_or(self, param1, param2):
        fetch_param_result = self._fetch_param_result

        def inside():
            result_param1 = fetch_param_result(param1)
            result_param2 = fetch_param_result(param2)

            return max(result_param1, result_param2)
        return inside

    def fuzzy_not(self, param1):
        # self.require_call_set_input_method()
        fetch_param_result = self._fetch_param_result

        def inside():
            result_param1 = fetch_param_result(param1)

            return 1 - result_param1
        return inside

    def set_input_data(self,input_data):
        self.input_data = input_data

    def _fetch_param_result(self,param1):
        if callable(param1):
            result_param1 = param1()
        elif type(param1) is str:
            field, member_f_name = param1.split()
            member_f = self.get_member_func(field,member_f_name)
            result_param1 = member_f.evaluate(self.input_data[field])
        else:
            result_param1 = param1
        return result_param1

    def mamdani(self, y: int, member_func: MemberFunct):
        evaluate = lambda x: min(member_func.evaluate(x), y)
        member_func._update_points(y)
        new_member_f = MemberFunct(member_func.field_name,member_func.memb_name,
                                   member_func.triangle_func, evaluate)
        new_member_f.points = member_func.points
        return new_member_f

    def larsen(self, y: int, member_func: MemberFunct):
        # member_func = self.get_member_func()
        evaluate = lambda x: y*member_func.evaluate(x)
        points = [(x, member_func.evaluate(x)) for x, _ in member_func.points]

        new_member_f = MemberFunct(member_func.field_name, member_func.memb_name,
                                   member_func.triangle_func, evaluate)
        new_member_f.points = points
        # este array ya está ordenado
        new_member_f.points.sort()
        return new_member_f

    def centroid(self, member_funct: MemberFunct):
        numer = quad(lambda x: x*member_funct.evaluate(x), 0, inf)
        denom = quad(member_funct.evaluate, 0, inf)
        return numer[0]/denom[0]

    @staticmethod
    def _get_max_point_neq_inf(points):
        points.sort()
        maxa = (0, 0)

        for i in range(len(points) - 1, 0, -1):
            maxa = points[i]
            if maxa[0] is not inf:
                break

        return maxa

    def lom(self, member_funct: MemberFunct):
        m = self._get_max_point_neq_inf(member_funct.points)
        a = [(i, j) for i, j in member_funct.points if j == m[1] and i != inf]
        mi = max(a)
        return mi[0]

    def som(self, member_funct: MemberFunct):
        m = self._get_max_point_neq_inf(member_funct.points)
        a = [(i, j) for i, j in member_funct.points if j == m[1] and i != inf]
        mi = min(a)
        return mi[0]

    def mom(self, member_funct: MemberFunct):
        m = self._get_max_point_neq_inf(member_funct.points)
        a = [i for i, j in member_funct.points if j == m[1] and i != inf]

        return sum(a)/len(a)

    def bisec(self, member_function: MemberFunct):

        ma = self._get_max_point_neq_inf(member_function.points)[0]
        epsilon = 1e-3
        if ma is inf:
            print('El área de la función de pertenencia es infinita')
            raise ValueError
        mi = min(member_function.points)[0]
        return self.binary_search(member_function.evaluate,ma,mi,mi,ma,ma/2,epsilon)

    def binary_search(self,member_function,max_value,min_value,lower_bound, uppper_bound,i,epsilon):
        area1 = quad(member_function,min_value,i)[0]
        area2 = quad(member_function,i,max_value)[0]

        r = area1 - area2
        if abs(r) < epsilon:
            return i
        # me muevo para la izq
        elif r > 0:
            # si

            return self.binary_search(
                    member_function,
                    max_value,
                    min_value,
                    lower_bound,
                    i,
                    (i + lower_bound) / 2,
                    epsilon)

        return self.binary_search(
                    member_function,
                    max_value,
                    min_value,
                    i,
                    uppper_bound,
                    (i + uppper_bound) / 2,
                    epsilon)
