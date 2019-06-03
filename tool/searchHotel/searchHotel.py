#! /usr/bin/python
# -*-coding: utf-8-*-

from selenium import webdriver
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait
import time
import datetime

from bs4 import BeautifulSoup
import codecs


class QunaSpider():
    def __init__(self):
        print("去哪儿网酒店信息")



    def get_hotel(self, dr, to_city, fromdate, todate):
        '''

        :param dr:
        :param to_city:
        :param fromdate:
        :param todate:
        :return:
        '''
        # ele_fromCity = dr.find_element_by_name("fromCity")
        #ele_fromDate = dr.find_element_by_id("fromDate")

        ele_toCity = dr.find_element_by_name("toCity")
        ele_toDate = dr.find_element_by_id("toDate")
        ele_search = dr.find_element_by_class_name("search-btn")

        ele_toCity.clear()
        ele_toCity.send_keys(to_city)
        ele_toCity.click()

        #ele_fromDate.clear()
        #ele_fromDate.send_keys(fromdate)

        ele_toDate.clear()
        ele_toDate.send_keys(todate)

        ele_search.click()

        page_num = 0

        while True:
            try:
                WebDriverWait(dr, 10).until(
                    EC.title_contains(unicode(to_city))
                )
            except Exception, e:
                print(e)
                break

            time.sleep(5)

            js = "window.scrollTo(0, document.body.scrollHeight);"
            dr.execute_script(js)
            time.sleep(5)

            html_const = dr.page_source
            soup = BeautifulSoup(html_const, 'html.parser', from_encoding='utf-8')
            infos = soup.find_all(class_="item_hotel_info")

            f = codecs.open(unicode(to_city) + unicode(fromdate) + u'.html', 'a+', 'utf-8')

            for info in infos:
                s1 = '=='*25 + str(page_num) + '=='*25
                print(s1)
                f.write(s1)
                content = info.get_text().replace(" ", "").replace("\t", "").strip()
                #print(content)
                for line in [ln for ln in content.splitlines() if ln.strip()]:
                    f.write(line)
                    f.write('\r\n')

            f.close()

            if page_num < 10:
                try:
                    next_page = WebDriverWait(dr, 10).until(
                        EC.visibility_of(dr.find_element_by_css_selector(".item.next"))
                    )
                    next_page.click()
                    page_num += 1
                    time.sleep(10)
                except Exception, e:
                    print(e)
                    break
            else:
                break


    def crawl(self, root_url, to_city):
        '''

        :param root_url:
        :param to_city:
        :return:
        '''
        today = datetime.date.today().strftime('%Y-%m-%d')
        tomorrow = datetime.date.today() + datetime.timedelta(days=1)
        tomorrow = tomorrow.strftime("%Y-%m-%d")

        dr = webdriver.Chrome(executable_path=r"/Users/software/chromedriver.exe")
        dr.set_page_load_timeout(50)
        dr.get(root_url)
        dr.maximize_window()    # 浏览器最大化显示
        dr.implicitly_wait(10)  # 控制间隔时间，等待浏览器反应
        self.get_hotel(dr, to_city, today, tomorrow)





if __name__ == '__main__':
    spider = QunaSpider()
    root_url = 'http://hotel.qunar.com/'
    spider.crawl(root_url, u'深圳')
