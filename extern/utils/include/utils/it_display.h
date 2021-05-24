/**
 * @file it_display.h
 * @brief Iteration Display Class
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 */

#ifndef IT_DISPLAY_H
#define IT_DISPLAY_H

#include <string>
#include <map>
#include <memory>
#include <fmt/format.h>


namespace dominiqs {

/**
 * Class for printing table-like iteration logs to the terminal
 */

class IterationDisplay
{
public:
	IterationDisplay() : headerInterval(100), iterationInterval(10), marked(false) {}
	~IterationDisplay() { clear(); }
	bool addColumn(const std::string& name, int priority=0, int width=10, int precision=2, bool visible=true, const std::string& d="-");
	void removeColumn(const std::string& name);
	void clear();
	void setVisible(const std::string& name, bool visible);
	void printHeader(std::ostream& out);
	bool needHeader(int k) const { return ((k % headerInterval) == 0); }
	void resetIteration();
	void markIteration() { marked = true; }
	bool needPrint(int k) const { return (marked || ((k % iterationInterval) == 0)); }
	template<typename T> void set(const std::string& name, const T& data)
	{
		bool found = false;
		for (auto& citr: columns)
		{
			ColumnPtr col = citr.second;
			if (col->name == name)
			{
				current[name] = col->formatValue(data);
				found = true;
				break;
			}
		}
		if (!found)  throw std::runtime_error(fmt::format("IterationDisplay: column {} does not exist", name));
	}
	void printIteration(std::ostream& out);
	int headerInterval;
	int iterationInterval;
protected:
	class Column
	{
	public:
		Column(const std::string& n, int p, int w, int prec, bool v, const std::string& d) : name(n), priority(p), width(w), precision(prec), visible(v), defValue(d) {}
		std::string name;
		int priority;
		int width;
		int precision;
		bool visible;
		std::string defValue;
		std::string formatHeader(const std::string& name) const
		{
			return fmt::format("{:>{}}", name, width);
		}
		template<typename T> std::string formatValue(const T& data) const
		{
			return fmt::format("{:>{}}", data, width);
		}
		std::string formatValue(const double& data) const
		{
			return fmt::format("{:>{}.{}f}", data, width, precision);
		}
	};
	using ColumnPtr = std::shared_ptr<Column>;
	using ColumnMap = std::map<int, ColumnPtr>;
	ColumnMap columns;
	using ItMap = std::map<std::string, std::string>;
	ItMap current;
	bool marked;
};


} // namespace dominiqs

#endif /* IT_DISPLAY_H */
